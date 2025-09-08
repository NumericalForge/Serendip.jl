# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

# Parallel loop with threads and the option of initializers for accumulators
macro withthreads_batch(ex)
    exbk = copy(ex)
    ex = Base.remove_linenums!(ex)

    if ex.head===:block # initializers and for loop
        # looking for initializers at the beginning: R=Int[]  or  R, C = Int[], Int[]
        initializers = Expr[]
        for arg in ex.args[1:end-1]
            if arg.head != :(=)
                throw(ArgumentError("@withthreads in a block requires initializers"))
            end
            if isa(arg.args[1], Symbol)
                push!(initializers, arg)
            elseif arg.args[1].head == :tuple
                vars = arg.args[1].args
                vals = arg.args[2].args
                for (var, val) in zip(vars,vals)
                    push!(initializers, :($var=$val))
                end
            else
                throw(ArgumentError("@withthreads bad initializer: $(arg)"))
            end

        end
        loop = ex.args[end]
        ex = Expr(:block)
        append!(ex.args, initializers)
        push!(ex.args, loop)

        # generating expressions and replacing inner variables
        initializers = ex.args[1:end-1]
        init_exp     = Expr(:block)
        reduce_exp   = Expr(:block)
        locals_exp   = Expr(:local)
        fetchvar     = Symbol("#_res")

        for (i,ini) in enumerate(initializers)
            var = ini.args[1]
            newvar = Symbol("#_$var")
            ex = replace(ex, var => newvar)

            push!(locals_exp.args, var)

            inilenvar = Symbol("#_$(var)_ini")
            push!(init_exp.args, ini )
            push!(init_exp.args, :( $inilenvar = length($var) ))
            push!(reduce_exp.args,
                :(
                    if $inilenvar==0
                        append!($var, $fetchvar[$i])
                    else
                        $var .+= $fetchvar[$i]
                    end
                )
            )
        end

        # generating inner expressions
        inner_initializers = ex.args[1:end-1]
        inner_return_exp = Expr(:return, :())
        inner_init_exp = Expr(:block)
        for ini in inner_initializers
            var = ini.args[1]
            push!(inner_return_exp.args[1].args, var)
            push!(inner_init_exp.args, ini)
        end

        # updating loop and body
        loop = ex.args[end]
        body = loop.args[2]
    else # only for loop
        loop = ex
    end

    if !(isa(loop, Expr) && loop.head === :for)
        throw(ArgumentError("@withthreads requires a `for` loop expression"))
    end
    itr  = loop.args[1].args[1]
    list = loop.args[1].args[2]

    return quote
        nt = Threads.nthreads()
        n  = length($(esc(list)))
        n<nt && (nt=n)
        if nt==1
            $(esc(exbk))
        else
            tasks = Vector{Task}(undef, nt)

            ccall(:jl_enter_threaded_region, Cvoid, ())

            # new task and schedule
            for tid in 1:nt
                function fun()
                    nmin = (tid-1)*div(n, nt) + 1
                    nmax = tid==nt ? n : tid*div(n, nt)
                    $(esc(inner_init_exp))
                    for $(esc(itr)) in $(esc(list))[nmin:nmax]
                        $(esc(body))
                    end
                    return $(esc(inner_return_exp))
                end

                t = Task(fun)
                t.sticky = true
                tasks[tid] = t
                ccall(:jl_set_task_tid, Cvoid, (Any, Cint), t, tid-1)
                schedule(t)
            end

            # gathering results
            $(esc(locals_exp))
            try
                $(esc(init_exp))
                for tid in 1:nt
                    $(esc(fetchvar)) = fetch(tasks[tid])
                    $(esc(reduce_exp))
                end
            finally
                ccall(:jl_exit_threaded_region, Cvoid, ())
            end
        end
    end
end


macro withthreads(ex)
    exbk = copy(ex)
    ex = Base.remove_linenums!(ex)

    if ex.head === :block
        # collect initializers e.g. R=Int[]  or  R, C = Int[], Int[]
        initializers = Expr[]
        for arg in ex.args[1:end-1]
            if arg.head != :(=)
                throw(ArgumentError("@withthreads in a block requires initializers"))
            end
            if isa(arg.args[1], Symbol)
                push!(initializers, arg)
            elseif arg.args[1].head == :tuple
                vars = arg.args[1].args
                vals = arg.args[2].args
                for (var, val) in zip(vars, vals)
                    push!(initializers, :($var = $val))
                end
            else
                throw(ArgumentError("@withthreads bad initializer: $(arg)"))
            end
        end
        loop = ex.args[end]
        ex = Expr(:block)
        append!(ex.args, initializers)
        push!(ex.args, loop)

        # prepare renaming and reductions
        initializers = ex.args[1:end-1]
        init_exp   = Expr(:block)
        reduce_exp = Expr(:block)
        locals_exp = Expr(:local)
        fetchvar   = Symbol("#_res")

        for (i, ini) in enumerate(initializers)
            var    = ini.args[1]
            newvar = Symbol("#_$var")
            ex     = replace(ex, var => newvar)
            push!(locals_exp.args, var)

            inilenvar = Symbol("#_$(var)_ini")
            push!(init_exp.args, ini)
            push!(init_exp.args, :( $inilenvar = length($var) ))
            push!(reduce_exp.args, :(
                if $inilenvar == 0
                    append!($var, $fetchvar[$i])
                else
                    $var .+= $fetchvar[$i]
                end
            ))
        end

        # inner returns and per-task init
        inner_initializers = ex.args[1:end-1]
        inner_return_exp = Expr(:return, :())
        inner_init_exp   = Expr(:block)
        for ini in inner_initializers
            var = ini.args[1]
            push!(inner_return_exp.args[1].args, var)
            push!(inner_init_exp.args, ini)
        end

        loop = ex.args[end]
        body = loop.args[2]
    else
        loop = ex
    end

    if !(isa(loop, Expr) && loop.head === :for)
        throw(ArgumentError("@withthreads requires a `for` loop expression"))
    end
    itr  = loop.args[1].args[1]
    list = loop.args[1].args[2]

    return quote
        nt = Threads.nthreads()
        n  = length($(esc(list)))
        nt = min(nt, max(n, 1))
        if nt == 1
            $(esc(exbk))
        else
            tasks = Task[]

            # allocate outer accumulators and capture their initial length
            $(esc(locals_exp))
            $(esc(init_exp))

            # dynamic chunking to handle uneven work
            # heuristic: ~8 chunks per thread
            chunk = max(1, cld(n, 8*nt))

            for s in 1:chunk:n
                e = min(s + chunk - 1, n)
                push!(tasks, Threads.@spawn begin
                    $(esc(inner_init_exp))
                    @inbounds for $(esc(itr)) in $(esc(list))[s:e]
                        $(esc(body))
                    end
                    $(esc(inner_return_exp))
                end)
            end

            for t in tasks
                $(esc(fetchvar)) = fetch(t)
                $(esc(reduce_exp))
            end
        end
    end
end
