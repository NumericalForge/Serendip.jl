# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

mutable struct Monitor{T}
    kind    ::Symbol
    exprs   ::Vector
    selector::Any
    target  ::Vector{T}
    stops   ::Vector
    filename::String
    helper_values::OrderedDict{Symbol, Number}
    values  ::OrderedDict{Union{Symbol,Expr}, Number}
    table   ::DataTable
end

_monitor_expr_name(expr::Symbol) = expr
_monitor_expr_name(expr::Expr) = expr.head == :(=) ? expr.args[1] : expr

function _monitor_expr_value(expr, state)
    if expr isa Expr && expr.head == :(=)
        var = expr.args[1]
        val = evaluate(expr.args[2]; state...)
        state[var] = val
        return val
    end
    return evaluate(expr; state...)
end

function _monitor_vars(exprs)
    vars = Symbol[]
    for expr in exprs
        append!(vars, getvars(expr))
    end
    return unique(vars)
end

function _monitor_minmax_base(var::Symbol)
    name = String(var)
    if endswith(name, "_min")
        return Symbol(replace(name, r"_min$" => "")), min
    elseif endswith(name, "_max")
        return Symbol(replace(name, r"_max$" => "")), max
    else
        return nothing, nothing
    end
end

function _monitor_minmax_values!(state, monitor::Monitor)
    vars = union(_monitor_vars(monitor.exprs), _monitor_vars(monitor.stops))
    vals = OrderedDict{Symbol, Number}()

    for var in vars
        base, fun = _monitor_minmax_base(var)
        fun === nothing && continue

        haskey(state, base) || error("Monitor: Value `$(base)` required by `$(var)` was not found. Available values are: $(join(keys(state), ", ")).")
        current = state[base]

        if !(var in keys(monitor.helper_values))
            monitor.helper_values[var] = current
        else
            monitor.helper_values[var] = fun(current, monitor.helper_values[var])
        end
        state[var] = monitor.helper_values[var]
        vals[var] = monitor.helper_values[var]
    end

    return vals
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
add_monitor(ana, :node, (x==0,y==0), :ux, "monitor_node_ux.dat")

# Monitor stress components at a group of integration points
add_monitor(ana, :ipgroup, z>1.0, [:sxx, :syy], "stress_ipgroup.dat")

# Monitor a node and stop analysis if ux > 0.01 at that node
add_monitor(ana, :node, [0.0, 0.0, 0.0], :ux; stop=:(ux > 0.01))
```
"""
function add_monitor(
    ana::Analysis,
    kind::Symbol,
    selector,
    expr::Union{Symbol,Expr,Tuple,Array,Nothing} = nothing,
    filename="";
    stop::Union{Expr,Tuple,Array,Nothing} = nothing,
    )

    @check kind in (:node, :ip, :nodegroup, :ipgroup, :nodalreduce) "add_monitor: kind must be one of :node, :ip, :nodegroup, :ipgroup, :nodalreduce"

    if filename != ""
        formats = (".dat", ".tab", ".table", ".json")
        _, format = splitext(filename)
        @check format in formats "Monitors must have one of the following extensions: $(join(formats, ", ", " and ")). Got $(repr(format))."
    end

    item_kind = kind in (:node, :nodegroup, :nodalreduce) ? :node : :ip
    target_type = item_kind == :node ? Node : Ip

    # get filename
    filename = getfullpath(ana.data.outdir, filename)

    exprs = if expr === nothing
        Any[]
    elseif expr isa Symbol
        [ expr ]
    elseif expr isa Expr 
        if expr.head == :tuple
            vec(expr.args)
        else
            [ expr ]
        end
    else        
        vec(collect(expr))
    end 

    isempty(exprs) && stop === nothing && error("add_monitor: provide at least one expression or one stop condition")

    all(typeof(ex) in (Symbol, Expr) for ex in exprs) || error("add_monitor: expr must be a Symbol or Expr or a collection of them")

    for (i,ex) in enumerate(exprs)
        ex = round_floats!(ex)
        ex = fix_expr_maximum_minimum!(ex)
        exprs[i] = ex
    end

    stop !== nothing && kind in (:nodegroup, :ipgroup) && error("Stop condition not supported for $kind monitors")
    stops = stop === nothing ? [] : stop isa Expr ? [stop] : vec(collect(stop))
    all(stop isa Expr for stop in stops) || error("add_monitor: stop must be a condition Expr or a collection of Exprs")

    for (i,expr) in enumerate(stops)
        ex = round_floats!(expr)
        ex = fix_expr_maximum_minimum!(ex)
        if ex isa Expr && ex.head == :(=)
            lhs = ex.args[1]
            rhs = ex.args[2]
            notify("add_monitor: Ignoring alias `$(lhs)` in stop condition")
            ex = rhs
        end
        stops[i] = ex
    end

    selector_str = selector isa String || selector isa Symbol ? repr(selector) : replace(string(selector), r"(?<!\,)\s+" => "")

    use_nearest = kind in (:node, :ip)
    target = select(ana.model, item_kind, selector; nearest=use_nearest, prefix="add_monitor")
    n = length(target)
    n == 0 && notify("add_monitor: No $(item_kind)s match: $selector_str")
    if kind in (:node, :ip)
        n >  1 && notify("add_monitor: Multiple $(item_kind)s match: $selector_str. Picking one")
        n >= 1 && (target = target[1:1])
    end

    mon = Monitor{target_type}(kind, exprs, selector, target, stops, filename, OrderedDict{Symbol, Number}(), OrderedDict{Union{Symbol,Expr}, Number}(), DataTable())
    update_monitor!(ana, mon)
    push!(ana.data.monitors, mon)
    return mon
end


function update_monitor!(ana::Analysis, monitor::Monitor)
    length(monitor.target) == 0 && return success()

    if monitor.kind in (:node, :ip)
        state = get_values(monitor.target[1])
        empty!(monitor.values)

        _monitor_minmax_values!(state, monitor)

        
        for expr in monitor.exprs
            if expr isa Symbol && !(expr in keys(state))
                error("Monitor: Value `$(expr)` not found at $(monitor.kind) $(monitor.target[1].id). Available values for monitoring are: $(join(keys(state), ", ")).")
            end
            monitor.values[_monitor_expr_name(expr)] = _monitor_expr_value(expr, state)
        end

        monitor.values[:stage] = ana.data.stage
        monitor.values[:T]     = ana.data.T
        ana.data.transient && (monitor.values[:t]=ana.data.t)
        push!(monitor.table, monitor.values)

        # eval stop expressions
        for expr in monitor.stops
            if evaluate(expr; state...)
                return failure("Analysis stopped by monitor condition: ($expr)")
            end
        end
    elseif monitor.kind in (:ipgroup, :nodegroup)
        empty!(monitor.values)
        for expr in monitor.exprs
            vals = []
            for item in monitor.target
                state = get_values(item)
                _monitor_minmax_values!(state, monitor)
                val = _monitor_expr_value(expr, state)
                push!(vals, val)
            end
            name = _monitor_expr_name(expr)
            if name isa Symbol && expr isa Symbol
                lo, hi = round.(extrema(vals), sigdigits=5)
                monitor.values[Symbol("$(name)_min")] = lo
                monitor.values[Symbol("$(name)_max")] = hi
            else
                val = any(vals) # TODO ?
                monitor.values[name] = val
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

        empty!(monitor.values)
        _monitor_minmax_values!(state, monitor)

        # eval expressions
        for expr in monitor.exprs
            monitor.values[_monitor_expr_name(expr)] = _monitor_expr_value(expr, state)
        end

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
        push!(labels, replace(string(k), " " => ""))
        push!(values, string(v))

    end
    return head, labels, values
end


function save(monitor::Monitor, filename::String; quiet=false)
    save(monitor.table, filename, quiet=quiet)
end
