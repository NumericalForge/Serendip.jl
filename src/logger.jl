# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

mutable struct Logger{T}
    kind    ::Symbol
    selector::Any
    target  ::Vector{T}
    filename::String
    table   ::DataTable
end


"""
    add_logger(ana::Analysis, kind::Symbol, selector::Any, filename::String = "")

Adds a logger to the analysis.

# Positional Arguments

- `ana`: The analysis object to add the logger to.
- `kind`: Type of logger to create. Can be one of `:node`, `:ip`, `:nodegroup`, `:ipgroup`, or `:nodalreduce`.
- `selector`: Filter expression to select the items to log. Can be a vector, a condition, or an expression.
- `filename`: Filename to save the log.
"""

"""
    add_logger(ana::Analysis, kind::Symbol, selector::Any, filename::String = "")

Adds a logger to an analysis, allowing you to record data during simulation.

# Arguments

- `ana::Analysis`:
    The analysis object to which the logger will be attached. The logger is stored in `ana.loggers`.

- `kind::Symbol`:
    Specifies the type of logger to create. Valid options are:
    - `:node`: Logs data at a single node.
    - `:ip`: Logs data at a single integration point.
    - `:nodegroup`: Logs data for a group of nodes.
    - `:ipgroup`: Logs data for a group of integration points.
    - `:nodalreduce`: Logs aggregated data (e.g., sum or average) across selected nodes.

- `selector::Any`:
    Defines how to select the items to log. Can be:
    - A **vector** `[x, y, z]` specifying coordinates (the nearest matching item will be used if no exact match is found).
    - A **logical expression** (e.g., `x > 0`).
    - A **predefined tag** or list of items.

- `filename::String` (optional):
    Name of the file where the log will be saved. If not provided, the logger will use the default output directory in `ana.data.outdir`.
    The file extension should be `.table`, `.table` or `.json`.

# Returns

- A `Logger` object of type `Logger{Node}` or `Logger{Ip}` depending on `kind`. This logger is also pushed into `ana.loggers`.

# Notes

- If `kind` is `:node` or `:ip` and the selector matches multiple items, only the **first match** will be logged.
- For `:nodegroup` and `:ipgroup`, the selected items are automatically sorted based on coordinates.
- If no items match the selector, a notification is displayed.
"""
function add_logger(
    ana::Analysis,
    kind::Symbol,
    selector::Any,
    filename::String = ""
)
    @check kind in (:node, :ip, :nodegroup, :ipgroup, :nodalreduce) "Unknown logger kind: $kind. Supported kinds are :node, :ip, :nodegroup, :ipgroup, :nodalreduce"
    
    # if filename != ""
    #     if kind in (:node, :ip, :nodalreduce) 
    #         formats = (".table", ".json")
    #     else
    #         formats = (".book", ".json")
    #     end
    #     _, format = splitext(filename)
    #     @check format in formats "Logger of kind $kind must have one of the following extensions: $(join(formats, ", ", " and ")). Got $(repr(format))."
    # end

    if filename != ""
        formats = (".tab", ".table", ".json")
        _, format = splitext(filename)
        @check format in formats "Loggers must have one of the following extensions: $(join(formats, ", ", " and ")). Got $(repr(format))."
    end
8
    item_name = kind in (:node, :nodegroup, :nodalreduce) ? :node : :ip
    target_type = item_name == :node ? Node : Ip

    filename = getfullpath(ana.data.outdir, filename)
    items = item_name == :node ? ana.model.nodes : select(ana.model, :ip, :all)

    if kind == :node
        name = "Node logger"
    elseif kind == :ip
        name = "Integration point logger"
    elseif kind == :nodalreduce
        name = "Nodal reduction logger"
    else 
        name = ""
    end

    if kind in (:node, :ip) && selector isa AbstractArray
        X = Vec3(selector)
        x, y, z = X
        target = select(ana.model, item_name, :(x==$x && y==$y && z==$z), nearest=false)
        n = length(target)
        if n==0
            notify("add_logger: No $kind found at $(selector). Picking the nearest at $X")
            target = [ nearest(items, X) ]
        else
            target = target[1:1] # take the first
        end

        log = Logger{target_type}(kind, selector, target, filename, DataTable(name=name))
        push!(ana.data.loggers, log)
        return log
    end

    target = select(ana.model, item_name, selector)
    n = length(target)
    n == 0 && notify("add_logger: No $(item_name)s found for selector: ", selector)
    if kind in (:node, :ip)
        n >  1 && notify("add_logger: More than one $item_name match selector: ", repr(selector), ". Picking the first one.")
        n >= 1 && (target = target[1:1])
    end

    if kind in (:nodegroup, :ipgroup)
        sort!(target, by=item->sum(item.coord))
    end

    log = Logger{target_type}(kind, selector, target, filename, DataTable(name=name))
    push!(ana.data.loggers, log)

    return log
end


function update_logger!(logger::Logger, ana::Analysis)
    length(logger.target) == 0 && return
    T = round(ana.data.stage - 1 + ana.data.T, digits=6)

    if logger.kind  in (:node, :ip)
        vals = get_values(logger.target[1])
        ana.data.transient && (vals[:t] = ana.data.t)
        push!(logger.table, vals)
    elseif logger.kind == :ipgroup
        # table = DataTable(name = "Ip group logger: T = $T")
        for ip in logger.target
            vals = get_values(ip)
            vals[:T] = T
            push!(logger.table, vals)
        end
        # push!(logger.book, table)
    elseif logger.kind == :nodegroup
        table = DataTable(name = "Node group logger: T = $T")

        # get data from recovered nodal values
        ids = [ node.id for node in logger.target ]
        for (key,V) in ana.model.node_data
            table[key] = V[ids]
        end

        # add coordinates
        for (i,key) in enumerate((:x, :y, :z))
            table[key] = [ item.coord[i] for item in logger.target ]
        end

        table[:T] = fill(T, length(logger.target))
        append!(logger.table, table)

        # push!(logger.book, table)
    elseif logger.kind == :nodalreduce
        tableU = DataTable()
        tableF = DataTable()
        for node in logger.target
            nvals = get_values(node)
            valsU = OrderedDict( dof.name => nvals[dof.name] for dof in node.dofs )
            valsF = OrderedDict( dof.natname => nvals[dof.natname] for dof in node.dofs )
            push!(tableF, valsF)
            push!(tableU, valsU)
        end

        valsU = OrderedDict( Symbol(key) => mean(tableU[key]) for key in keys(tableU) ) # gets the average of essential values
        valsF = OrderedDict( Symbol(key) => sum(tableF[key])  for key in keys(tableF) ) # gets the sum for each component
        vals  = merge(valsU, valsF)
        vals[:out] = ana.data.out
        ana.data.transient && (vals[:t] = ana.data.t)

        push!(logger.table, vals)
    end

end


function save(logger::Logger, filename::String; quiet=false)
    save(logger.table, filename, quiet=quiet)
    # if kind in (:node, :ip, :nodalreduce)
    # else
    #     save(logger.book, filename, quiet=quiet)
    # end
end
