# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

"""
    @withthreads block

A macro for parallelizing `for` loops with support for thread-local accumulators 
and a fail-fast mechanism.

# Arguments
The macro expects a `begin ... end` block containing:
1. **Initializers**: One or more assignments (e.g., `R = Int[]`) defining variables 
   to be accumulated.
2. **For Loop**: A standard Julia `for` loop.

# Behavior
- **Parallelization**: Splits the loop iterations across available `Threads.nthreads()`.
- **Accumulation**: Each thread works on a private copy of the initializers. 
  - If the initial variable is empty (length 0), results are combined using `append!`.
  - Otherwise, results are combined using in-place addition `.+=`.
- **Fail-Fast**: 
  - Monitors a shared `Threads.Atomic{Bool}` stop signal.
  - If any task encounters an exception or executes a `break` statement, all other 
    tasks will terminate at the start of their next iteration.
- **Exception Handling**: Re-throws any exception encountered during execution to 
  the main thread.

# Example: Finite Element Assembly
```
@withthreads begin
    R, C, V = Int[], Int[], Float64[]
    for elem in elements
        ke, rmap, cmap = compute_stiffness(elem)
        append!(R, rmap); append!(C, cmap); append!(V, ke)
    end
end
```
"""
macro withthreads_new_slower(ex)
    exbk = copy(ex)
    ex = Base.remove_linenums!(ex)

    # Recursive helper to replace symbols and expressions
    function walk(f, x)
        if x isa Expr
            return f(Expr(x.head, map(arg -> walk(f, arg), x.args)...))
        else
            return f(x)
        end
    end

    if ex.head === :block # the macro expectas a block with initializers and a for loop
        initializers = Expr[]
        for arg in ex.args[1:end-1]
            if arg.head != :(=)
                throw(ArgumentError("@withthreads in a block requires initializers"))
            end
            if isa(arg.args[1], Symbol)
                push!(initializers, arg)
            elseif arg.args[1].head == :tuple
                vars, vals = arg.args[1].args, arg.args[2].args
                for (var, val) in zip(vars, vals)
                    push!(initializers, :($var = $val))
                end
            end
        end
        loop = ex.args[end] # the loop code
        
        init_exp         = Expr(:block)
        reduce_exp       = Expr(:block)
        locals_exp       = Expr(:local)
        fetchvar         = Symbol("#_res")
        stop_signal      = Symbol("#_stop_signal")
        inner_init_exp   = Expr(:block)
        inner_return_exp = Expr(:return, Expr(:tuple))

        for (i, ini) in enumerate(initializers)
            var    = ini.args[1]
            newvar = Symbol("#_$var")
            
            # Replace variable with thread-local version
            loop = walk(x -> x === var ? newvar : x, loop)

            push!(locals_exp.args, var)
            ini_len_var = Symbol("#_$(var)_ini")
            push!(init_exp.args, ini)
            push!(init_exp.args, :($ini_len_var = length($var)))
            push!(reduce_exp.args, 
                :(if $ini_len_var == 0
                    append!($var, $fetchvar[$i])
                  else
                    $var .+= $fetchvar[$i]
                  end)
            ) # accumulation/reduction
            push!(inner_init_exp.args, :($newvar = $(ini.args[2])))
            push!(inner_return_exp.args[1].args, newvar)
        end
        
        body = loop.args[2] # loop body
        itr  = loop.args[1].args[1] # iterator
        list = loop.args[1].args[2] # list of elements
    else
        throw(ArgumentError("@withthreads requires a block with initializers"))
    end

    # Replace 'break' with 'signal + break'
    body = walk(body) do x
        if x isa Expr && x.head === :break
            return quote
                $stop_signal[] = true
                break
            end
        end
        return x
    end

    return quote
        nt = Threads.nthreads()
        n  = length($(esc(list)))
        n < nt && (nt = n)
        
        if nt == 1
            $(esc(exbk))
        else
            tasks = Vector{Task}(undef, nt)
            $(esc(stop_signal)) = Threads.Atomic{Bool}(false)

            for tid in 1:nt
                function fun()
                    nmin = (tid - 1) * div(n, nt) + 1
                    nmax = tid == nt ? n : tid * div(n, nt)
                    $(esc(inner_init_exp))
                    
                    for $(esc(itr)) in $(esc(list))[nmin:nmax]
                        if $(esc(stop_signal))[]
                            break
                        end
                        try
                            $(esc(body))
                        catch e
                            $(esc(stop_signal))[] = true
                            rethrow(e)
                        end
                    end
                    return $(esc(inner_return_exp))
                end

                t = Task(fun)
                # t.sticky = true
                tasks[tid] = t
                schedule(t)
            end

            $(esc(locals_exp))
            try
                $(esc(init_exp))
                for tid in 1:nt
                    $(esc(fetchvar)) = fetch(tasks[tid])
                    $(esc(reduce_exp))
                end
            catch err
                $(esc(stop_signal))[] = true
                rethrow(err)
            end
        end
    end
end


"""
    @withthreads block

A macro for parallelizing `for` loops with support for thread-local accumulators 
and a fail-fast mechanism.

# Arguments
The macro expects a `begin ... end` block containing:
1. **Initializers**: One or more assignments (e.g., `R = Int[]`) defining variables 
   to be accumulated.
2. **For Loop**: A standard Julia `for` loop.

# Behavior
- **Parallelization**: Splits the loop iterations across available `Threads.nthreads()`.
- **Accumulation**: Each thread works on a private copy of the initializers. 
  - If the initial variable is empty (length 0), results are combined using `append!`.
  - Otherwise, results are combined using in-place addition `.+=`.
- **Fail-Fast**: 
  - Monitors a shared `Threads.Atomic{Bool}` stop signal.
  - If any task encounters an exception or executes a `break` statement, all other 
    tasks will terminate at the start of their next iteration.
- **Exception Handling**: Re-throws any exception encountered during execution to 
  the main thread.

# Example: Finite Element Assembly
```
@withthreads begin
    R, C, V = Int[], Int[], Float64[]
    for elem in elements
        ke, rmap, cmap = compute_stiffness(elem)
        # If an error occurs here and you 'break', all threads stop.
        if is_invalid(ke)
            break
        end
        append!(R, rmap); append!(C, cmap); append!(V, ke)
    end
end
```
"""
macro withthreads(ex)
    exbk = copy(ex)
    ex = Base.remove_linenums!(ex)

    if ex.head===:block # initializers and for loop
        # looking for initializers: R=Int[]  or  R, C = Int[], Int[]
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

        # generating expressions and replaceing inner variables
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

            # new task and schedule
            for tid in 1:nt
                nmin = (tid-1) * div(n, nt) + 1
                nmax = (tid == nt) ? n : tid * div(n, nt)

                tasks[tid] = Threads.@spawn begin
                    $(esc(inner_init_exp))
                    for k in nmin:nmax
                        $(esc(itr)) = $(esc(list))[k]
                        $(esc(body))
                    end
                    $(esc(inner_return_exp))
                end
            end

            # gathering results
            $(esc(locals_exp))
            try
                $(esc(init_exp))
                for tid in 1:nt
                    $(esc(fetchvar)) = fetch(tasks[tid])
                    $(esc(reduce_exp))
                end
            catch err
                rethrow(err)
            end
        end
    end
end





# Parallel loop with threads with option for accumulators 
# For some reason this macro accumulates greater errors during solving
macro withthreads_x(ex)
    exbk = copy(ex)
    ex   = Base.remove_linenums!(ex)

    # ❱❱❱ Parse block for initializers and for loop
    initializers = Expr[]
    loop = nothing
    if ex isa Expr && ex.head === :block
        isempty(ex.args) && throw(ArgumentError("@withthreads: block requires initializers and a for loop"))
        for arg in ex.args[1:end-1]
            (arg isa Expr && arg.head === :(=)) || throw(ArgumentError("@withthreads: block requires initializer assignments"))
            push!(initializers, arg)
        end
        loop = ex.args[end]
    else
        loop = ex
    end

    (loop isa Expr && loop.head === :for) || throw(ArgumentError("@withthreads: a `for` loop expression is required"))

    # ❱❱❱ Unpack the for loop
    itr  = loop.args[1].args[1]
    list = loop.args[1].args[2]
    body = loop.args[2]

    # ❱❱❱ Collect accumulators and initial values
    acc_vars = Symbol[]
    acc_vals = Any[]
    for ini in initializers
        lhs, rhs = ini.args
        if lhs isa Symbol
            push!(acc_vars, lhs)
            push!(acc_vals, rhs)
        elseif lhs isa Expr && lhs.head === :tuple
            vars = lhs.args
            vals = rhs.args
            length(vars) == length(vals) || throw(ArgumentError("@withthreads: tuple init arity mismatch"))
            for (v, val) in zip(vars, vals)
                v isa Symbol || throw(ArgumentError("@withthreads bad tuple initializer"))
                push!(acc_vars, v); push!(acc_vals, val)
            end
        else
            throw(ArgumentError("@withthreads: bad initializer"))
        end
    end

    # rename map for per-task accumulators: R -> __wt_R
    rename_map = Dict{Symbol,Symbol}(v => Symbol(:_wt_, v) for v in acc_vars)


    # ❱❱❱ Helper functions to detect and rewrite `break` statements
    function has_break(x)
        x isa Expr || return false
        x.head === :break && return true
        any(has_break, x.args)
    end
    hb = has_break(body)


    function rw(x)
        if x isa Expr
            if x.head === :break
                return :(_break_reached[] = true; break)
            else
                return Expr(x.head, map(rw, x.args)...)
            end
        elseif x isa Symbol
            return haskey(rename_map, x) ? rename_map[x] : esc(x)
        else
            return x
        end
    end
    body_rw = rw(body)


    # ❱❱❱ Build blocks

    # outer accumulators
    outer_init = Expr(:block, (:( $(esc(acc_vars[i])) = $(esc(acc_vals[i])) ) for i in eachindex(acc_vars))...)

    # initial lengths for reduction decision
    length_vec_sym = :_length_vec
    len_block = Expr(:block, :(local $length_vec_sym = Int[]),
                     (:( push!($length_vec_sym, length($(esc(acc_vars[i])))) ) for i in eachindex(acc_vars))...)

    # per-task locals (shadowed names)
    inner_init = Expr(:block, (:( $(rename_map[v]) = $(esc(acc_vals[i])) ) for (i,v) in pairs(acc_vars))...)

    # return tuple of per-task accumulators
    ret_tuple = Expr(:tuple, (rename_map[v] for v in acc_vars)...)

    # reduction block
    local_results_sym = :_local_results
    reduce_block = Expr(:block,
        (quote
            if $(length_vec_sym)[$i] == 0 # append 
                append!($(esc(acc_vars[i])), $local_results_sym[$i])
            else # accumulate
                $(esc(acc_vars[i])) .+= $local_results_sym[$i]
            end
         end for i in eachindex(acc_vars))...)

    # early-cancel guard
    guard_stmt = hb ? :( if _break_reached[]; break; end ) : :( nothing )

    # ❱❱❱ Final assembly

    return quote
        nt = Threads.nthreads()   # number of threads
        n  = length($(esc(list))) # number of elements
        nt = min(nt, max(n, 1))   # update number of threads 
        if nt == 1
            $(esc(exbk))
        else
            tasks = Task[]

            # outer accumulators and their base lengths
            $(outer_init)
            $(len_block)

            # shared cancellation flag
            _break_reached = Threads.Atomic{Bool}(false)

            # chunking heuristic
            max_chunks = nt*6
            chunk = max(1, cld(n, max_chunks))

            for s in 1:chunk:n
                e = min(s + chunk - 1, n)
                push!(tasks, Threads.@spawn begin
                    $(inner_init)
                    for $(esc(itr)) in $(esc(list))[s:e]
                        $(guard_stmt)
                        $(body_rw)
                    end
                    $(ret_tuple)
                end)
            end

            # for _t in tasks
            for i in eachindex(tasks)
                $local_results_sym = fetch(tasks[i])
                $(reduce_block)
            end
        end
    end
end
