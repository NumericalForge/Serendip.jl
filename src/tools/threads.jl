# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

# Parallel loop with threads with option for accumulators
macro withthreads(ex)
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
