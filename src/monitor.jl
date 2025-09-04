# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

mutable struct Monitor{T}
    kind    ::Symbol
    exprs   ::Vector
    selector::Any
    target  ::Vector{T}
    stops   ::Vector
    filename::String
    values  ::OrderedDict{Symbol, Float64}
    table   ::DataTable
end

"""
    add_monitor(
        ana::Analysis,
        kind::Symbol,
        selector::Any,
        expr::Union{Symbol, Expr, Tuple, Array},
        filename::String = "";
        stop::Union{Expr,Tuple,Array,Nothing} = nothing
    )

Adds a monitor to the analysis for observing quantities during the simulation.
Monitors are used to evaluate expressions at selected nodes or integration points and optionally trigger stopping conditions.

# Positional Arguments

- `ana::Analysis`:
    The analysis object to which the monitor will be attached. The monitor is stored in `ana.monitors`.

- `kind::Symbol`:
    Specifies the type of monitor. Valid options are:
    - `:node`: Monitor a single node.
    - `:ip`: Monitor a single integration point.
    - `:nodegroup`: Monitor a group of nodes.
    - `:ipgroup`: Monitor a group of integration points.
    - `:nodalreduce`: Monitor reduced quantities across selected nodes.

- `selector`:
    Defines how to select the items to monitor. Can be:
    - A **vector** `[x, y, z]` specifying coordinates (the nearest matching item will be used if no exact match is found).
    - A **logical expression** (e.g., `x > 0 && y < 1`).
    - A predefined tag or list of items.

- `expr`:
    One or more expressions to monitor. Can be a single `Symbol` (e.g., `:ux`), or a collection (e.g., `[:ux, :uy]`, a tuple, or an array).

- `filename::String` (optional):
    Name of the file where the monitor results will be saved. If a path is not
    provided, the file is saved in the current directory.
    The file extension should be `.table`, `.table` or `.json`.


# Keyword Arguments

- `stop` (optional):
    Defines stop conditions for the monitor. Can be:
    - An condition (e.g., `:(ux > 0.1)`).
    - A collection of conditions.
    **Note:** Stop conditions are only supported for `:node` and `:ip` kinds.

# Returns

- A `Monitor` object of type `Monitor{Node}` or `Monitor{Ip}`, depending on `kind`. The monitor is also pushed into `ana.monitors`.

# Notes

- If `kind` is `:node` or `:ip` and the selector matches multiple items, only the **first match** will be monitored.
- For `:nodegroup` and `:ipgroup`, the selected items are sorted automatically by coordinate sum.
- If no items match the selector, the nearest item is chosen and a notification is displayed.

# Example
```julia
# Monitor displacement ux at a node at x=0 and y==0
add_monitor(ana, :node, (x==0,y==0), :ux, "monitor_node_ux.table")

# Monitor stress components at a group of integration points
add_monitor(ana, :ipgroup, z>1.0, [:sxx, :syy], "stress_ipgroup.table")

# Monitor a node and stop analysis if ux > 0.01 at that node
add_monitor(ana, :node, [0.0, 0.0, 0.0], :ux; stop=:(ux > 0.01))
```
"""
function add_monitor(
    ana::Analysis,
    kind::Symbol,
    selector,
    expr::Union{Symbol,Expr,Tuple,Array},
    filename="";
    stop::Union{Expr,Tuple,Array,Nothing} = nothing,
    )

    @check kind in (:node, :ip, :nodegroup, :ipgroup, :nodalreduce) "add_monitor: kind must be one of :node, :ip, :nodegroup, :ipgroup, :nodalreduce"

    if filename != ""
        formats = (".tab", ".table", ".json")
        _, format = splitext(filename)
        @check format in formats "Monitors must have one of the following extensions: $(join(formats, ", ", " and ")). Got $(repr(format))."
    end

    item_kind = kind in (:node, :nodegroup, :nodalreduce) ? :node : :ip
    target_type = item_kind == :node ? Node : Ip

    # get filename
    filename = getfullpath(ana.data.outdir, filename)

    exprs = typeof(expr) in (Symbol, Expr) ? [expr] : vec(collect(expr))
    all(typeof(expr) in (Symbol, Expr) for ex in exprs) || error("add_monitor: expr must be a Symbol or Expr or a collection of them")

    for (i,ex) in enumerate(exprs)
        ex = round_floats!(ex)
        ex = fix_expr_maximum_minimum!(ex)
        exprs[i] = ex
    end

    stop !== nothing && kind in (:nodegroup, :ipgroup) && error("Stop condition not supported for $kind monitors")
    stops = stop === nothing ? [] : stop isa Expr ? [stop] : vec(collect(stop))
    all(stop isa Expr for stop in stops) || error("add_monitor: stop must be a condition Expr or a collection of Exprs")

    for (i,expr) in enumerate(stops)
        ex = round_floats!(ex)
        ex = fix_expr_maximum_minimum!(ex)
        stops[i] = ex
    end

    items = target_type == :node ? ana.model.nodes : select(ana.model, :ip, :all)

    # setup if selector is an array
    if kind in (:node, :ip) && selector isa AbstractArray
        X = Vec3(selector)
        x, y, z = X
        target = select(ana.model, item_kind, :(x==$x && y==$y && z==$z), nearest=false)

        n = length(target)
        if n==0
            target = [ nearest(items, X) ]
            X = target[1].coord
            notify("add_monitor: No $kind found at $(selector). Picking the nearest at $X")
        else
            target = target[1:1] # take the first
        end

        mon = Monitor{target_type}(kind, exprs, selector, target, stops, filename, OrderedDict{Symbol, Vector{Float64}}(), DataTable())
        update_monitor!(ana, mon)
        push!(ana.data.monitors, mon)
        return mon
    end

    target = select(ana.model, item_kind, selector)
    n = length(target)
    n == 0 && notify("setup_monitor: No $(item_kind)s found for selector: ", selector)
    if kind in (:node, :ip)
        n >  1 && notify("setup_monitor: More than one $item_kind match selector: ", selector, ". Picking the first one.")
        n >= 1 && (target = target[1:1])
    end

    mon = Monitor{target_type}(kind, exprs, selector, target, stops, filename, OrderedDict{Symbol, Vector{Float64}}(), DataTable())
    update_monitor!(ana, mon)
    push!(ana.data.monitors, mon)
    return mon
end


function update_monitor!(ana::Analysis, monitor::Monitor)
    length(monitor.target) == 0 && return success()

    if monitor.kind in (:node, :ip)
        state = get_values(monitor.target[1])

        # update state with _max and _min
        vars = union(getvars(monitor.exprs)..., getvars(monitor.stops))
        minmax_vals = OrderedDict()
        for var in vars
            contains(string(var), r"_max|_min") || continue
            keystr, optstr = split(string(var), '_')
            key = Symbol(keystr)

            if size(monitor.table,1) == 0
                state[var] = state[key]
            else
                fun = optstr=="min" ? min : max
                state[var] = fun(state[key], monitor.table[var][end])
            end
            minmax_vals[var] = state[var]
        end

        for expr in monitor.exprs
            monitor.values[expr] = evaluate(expr; state...)
        end

        merge!(monitor.values, minmax_vals)
        monitor.values[:stage] = ana.data.stage
        monitor.values[:T]     = ana.data.T
        ana.data.transient && (monitor.values[:t]=ana.data.t)
        push!(monitor.table, monitor.values)

        # eval stop expressions
        for expr in monitor.stops
            if evaluate(expr; state...)
                return failure("Stop condition reached: ($expr)")
            end
        end
    elseif monitor.kind in (:ipgroup, :nodegroup)
        for expr in monitor.exprs
            vals = []
            for item in monitor.target
                state = get_values(item)
                val = evaluate(expr; state...)
                push!(vals, val)
            end
            if expr isa Symbol
                lo, hi = round.(extrema(vals), sigdigits=5)
                monitor.values[Symbol(:(min($expr)))] = lo
                monitor.values[Symbol(:(max($expr)))] = hi
            else
                val = any(vals) # TODO ?
                monitor.values[expr] = val
            end
        end

        monitor.values[:stage] = ana.data.stage
        monitor.values[:T]     = ana.data.T
        ana.data.transient && (monitor.values[:t]=ana.data.t)

        push!(monitor.table, monitor.values)
    elseif monitor.kind == :nodalreduce
        tableU = DataTable()
        tableF = DataTable()
        for node in monitor.target
            nvals = get_values(node)
            valsU  = OrderedDict( dof.name => nvals[dof.name] for dof in node.dofs )
            valsF  = OrderedDict( dof.natname => nvals[dof.natname] for dof in node.dofs )
            push!(tableF, valsF)
            push!(tableU, valsU)
        end

        valsU = OrderedDict( Symbol(key) => mean(tableU[key]) for key in keys(tableU) ) # gets the average of essential values
        valsF = OrderedDict( Symbol(key) => sum(tableF[key])  for key in keys(tableF) ) # gets the sum for each component
        state = merge(valsU, valsF)

        # update state with definitinos
        for ex in monitor.exprs
            if ex isa Expr && ex.head == :(=)
                var = ex.args[1]
                state[var] = evaluate(ex.args[2]; state...)
            end
        end

        # update state with _max and _min
        vars = union(getvars(monitor.exprs), getvars(monitor.stops))
        minmax_vals = OrderedDict()
        for var in vars
            contains(string(var), r"_max|_min") || continue
            keystr, optstr = split(string(var), '_')
            key = Symbol(keystr)

            if size(monitor.table)[1] == 0
                state[var] = state[key]
            else
                fun = optstr=="min" ? min : max
                state[var] = fun(state[key], monitor.table[var][end])
            end
            minmax_vals[var] = state[var]
        end

        # eval expressions
        for expr in monitor.exprs
            if expr isa Expr && expr.head == :(=)
                expr = expr.args[1]
            end
            monitor.values[expr] = evaluate(expr; state...)
        end

        merge!(monitor.values, minmax_vals)
        monitor.values[:stage] = ana.data.stage
        monitor.values[:T]     = ana.data.T
        ana.data.transient && (monitor.values[:t]=ana.data.t)
        push!(monitor.table, monitor.values)

        # eval stop expressions
        for expr in monitor.stops
            if evaluate(expr; state...)
                return failure("Stop condition at NodeSumMonitor ($expr)")
            end
        end
    end

    return success()

end


function output(monitor::Monitor)
    selector = monitor.selector

    if selector isa Expr || selector isa Tuple
        loc = replace( string(round_floats!(getexpr(selector))), " " => "")
    elseif selector isa Integer
        loc = "$item_name $selector"
    else
        loc = repr(selector)
    end

    head = string(monitor.kind) * " at " * loc

    labels = String[]
    values = String[]
    
    for (k,v) in monitor.values
        k in (:stage, :T, :t) && continue
        if v isa Real && !(v isa Bool)
            v = round(v, sigdigits=5)
        end
        push!(labels, string(k))
        push!(values, string(v))

    end
    return head, labels, values
end


function save(monitor::Monitor, filename::String; quiet=false)
    save(monitor.table, filename, quiet=quiet)
end
