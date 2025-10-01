# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export DataTable, save, split_by, randtable
export compress, resize, filter, cut!, clamp!, denoise!
export getheader


# DataTable object
const KeyType = Union{Symbol,AbstractString}

# """
#     DataTable(header::AbstractArray, columns::Vector{<:AbstractVector}=Vector{Any}[]; name="")
#     DataTable(; pairs...)
#     DataTable(header::AbstractArray, matrix::Array{T,2} where T; name="")

# Construct a tabular data container with named columns.

# # Arguments
# - `header`: Array of column names (strings or convertible).
# - `columns`: Vector of column vectors. Defaults to empty vectors of matching length.
# - `name`: Optional table name (string).
# - `pairs...`: Keyword arguments mapping names to column vectors.
# - `matrix`: 2D array whose columns are converted into vectors.

# # Behavior
# - Enforces unique column names.
# - Builds an index mapping names → column positions.
# - Columns can be accessed as fields (`table.colname`) or by key.

# # Examples
# ```julia
# DataTable(["x","y"], [ [1,2], [3,4] ])
# DataTable(x=[1,2,3], y=[4,5,6])
# DataTable(["a","b"], [1 2; 3 4; 5 6])
# ```
# """
"""
    DataTable(header, columns=Vector[]; name="")
    DataTable(; pairs...)
    DataTable(header, matrix; name="")
    DataTable(filename, delim='\\t')

Construct a tabular container with named columns.

# Arguments
- `header::AbstractArray`: column names (strings or symbols).
- `columns::Vector{<:AbstractVector}`: one vector per column. All columns must have the same length.
- `name::AbstractString`: optional table name.
- `pairs...`: keyword form `colName = vector`. A special key `name` sets the table name.
- `matrix::AbstractMatrix`: a 2D array whose columns become table columns.
- `filename::AbstractString`: path to a delimited text file (`.dat`, `.table`, `.json`).
- `delim::Char`: column delimiter when reading from file (default `'\\t'`).

# Examples
```julia
DataTable(["x","y"], [[1,2,3], [10.0, 20.0, 30.0]]; name="XY")
DataTable(x=[1,2,3], y=[10.0,20.0,30.0], name="XY")
DataTable(["a","b"], [1 2; 3 4; 5 6]; name="M")
DataTable("data.table")                # reads header + data from file
```
"""
mutable struct DataTable
    columns::Vector{AbstractVector}
    colmap ::OrderedDict{String,Int} # Data index
    header ::Vector{String}
    name   ::String

    function DataTable(header::AbstractArray, columns::Vector{<:AbstractVector}=Vector{Any}[]; name::String="")
        @check length(header) == length(unique(header)) "DataTable: Header contains repeated keys"
        if length(columns)==0
            columns = AbstractVector[ [] for i in 1:length(header) ]
        end

        colmap = OrderedDict{String,Int}( string(key)=>i for (i,key) in enumerate(header) )

        return new(columns, colmap, string.(header), name)
    end

    function DataTable(; pairs...)
        columns = AbstractVector[]
        header  = KeyType[]
        name    = ""
        for (key, val) in pairs
            if key==:name && typeof(val) <: AbstractString
                name = string(val)
                continue
            end
            @check typeof(val) <: AbstractVector "DataTable: value for key $(repr(key)) is not a vector"
            push!(header, string(key))
            push!(columns, copy(val))
        end

        colmap = OrderedDict{String,Int}( string(key)=>i for (i,key) in enumerate(header) )

        return new(columns, colmap, string.(header), name)
    end

end

getcolumns(table::DataTable) = getfield(table, :columns)
getheader(table::DataTable)  = getfield(table, :header)
getcolidx(table::DataTable)  = getfield(table, :colmap)
getname(table::DataTable)    = getfield(table, :name)


function DataTable(header::AbstractArray, matrix::Array{T,2} where T; name::String="")
    nkeys = length(header)
    ncols = size(matrix,2)
    nkeys != ncols && error("DataTable: header and number of data columns do not match")
    types = [ typeof(matrix[1,i]) for i in 1:ncols ]
    columns = AbstractVector[ convert(Array{types[i],1}, matrix[:,i]) for i in 1:ncols ]
    return DataTable(header, columns, name=name)
end


function Base.size(table::DataTable)
    columns = getcolumns(table)
    ncols   = length(columns)
    nrows   = ncols==0 ? 0 : length(columns[1])
    return (nrows, ncols)
end


function Base.size(table::DataTable, dim::Int)
    columns = getcolumns(table)

    ncols   = length(columns)
    dim==2 && return ncols

    nrows = ncols==0 ? 0 : length(columns[1])
    dim==1 && return nrows

    error("size: Table dimension value should be 1 or 2. Got $dim")
end


function Base.getproperty(table::DataTable, name::Symbol)
    key   = string(name)
    colmap = getcolidx(table)
    if haskey(colmap, key)
        return getindex(table, key)
    end

    return getfield(table, name)
end


function Base.setproperty!(table::DataTable, sym::Symbol, val)
    key = string(sym)
    setindex!(table, val, key)
end


function Base.keys(table::DataTable)
    return keys(getcolidx(table))
end


function Base.push!(table::DataTable, row::Array{T,1} where T)
    nr, nc = size(table)
    @assert nc==length(row)
    columns = getcolumns(table)

    for (i,v) in enumerate(row)
        push!(columns[i], v)
    end
end


function Base.push!(table::DataTable, row::Union{NTuple, AbstractDict, Vector{<:Pair}, NTuple{N, <:Pair} where N})
    columns = getcolumns(table)
    colmap  = getcolidx(table)
    header  = getheader(table)
    nr, nc = size(table)

    if nc==0
        for pair in row
            key = string(pair.first)
            push!(header, key)
            push!(columns, [pair.second])
            colmap[key] = length(header)
        end
    else
        for pair in row
            key = string(pair.first)
            idx = get(colmap, key, 0)
            if idx==0
                # add new column
                push!(header, key)
                newcol = zeros(nr)
                push!(newcol, pair.second)
                push!(columns, newcol)
                colmap[key] = length(columns)
            else
                push!(columns[idx], pair.second)
            end
        end

        # add zero for missing values
        for col in columns
            if length(col)==nr
                if eltype(col)==String
                    push!(col, "")
                else
                    push!(col, 0.0)
                end
            end
        end
    end

    return row
end

function Base.append!(table::DataTable, table2::DataTable)
    nr1, nc1 = size(table)
    nr2, nc2 = size(table2)

    columns = getcolumns(table)
    colmap  = getcolidx(table)
    header  = getheader(table)

    for key in keys(table2)
        col2 = table2[key]
        idx1 = get(colmap, key, 0)
        if idx1==0
            # add new column
            push!(header, key)
            zr = eltype(col2)==String ? "" : zero(eltype(col2))
            newcol = fill(zr, nr1)
            append!(newcol, col2)
            push!(columns, newcol)
            colmap[key] = length(columns)
        else
            append!(table.columns[idx1], col2)
        end
    end

    # add zero for missing values
    nr = size(table,1)
    for col in table.columns
        if length(col)<nr
            zr = eltype(col)==String ? "" : zero(eltype(col))
            append!(col, fill(zr, nr-length(col)))
        end
    end

    return table
end


function Base.setindex!(table::DataTable, column::AbstractVector, key::KeyType)
    columns = getcolumns(table)
    colmap  = getcolidx(table)
    header  = getheader(table)

    if length(columns)>0
        n1 = length(columns[1])
        n2 = length(column)
        n1>0 && n1!=n2 && error("setindex! : vector length ($n2) for data ($key) is incompatible with DataTable rows ($n1)")
    end

    key = string(key)
    if haskey(colmap, key)
        idx = colmap[key]
        columns[idx] = column
    else
        push!(columns, column)
        push!(header, key)
        colmap[key] = length(columns)
    end
    return column
end


function Base.getindex(table::DataTable, key::KeyType)
    key    = string(key)
    colmap = getcolidx(table)
    haskey(colmap, key) || error("getindex: key $(repr(key)) not found in DataTable. Available keys are $(keys(table))")

    idx  = colmap[key]
    return getcolumns(table)[idx]
end


function Base.haskey(table::DataTable, key::KeyType)
    colmap = getcolidx(table)
    return haskey(colmap, string(key))
end


function Base.getindex(table::DataTable, keys::Array{<:KeyType,1})
    columns = [ table[string(key)] for key in keys ]
    return DataTable(keys, columns)
end


function Base.getindex(table::DataTable, rowidx::Int)
    columns = getcolumns(table)
    return [ col[rowidx] for col in columns ]
end


function Base.getindex(table::DataTable, idxs::Union{Colon,UnitRange{Int},Array{Int,1},BitArray{1}})
    columns = getcolumns(table)
    cols = [ columns[i][idxs] for i in 1:length(columns) ]

    return DataTable(getheader(table), cols)
end


function Base.getindex(table::DataTable, rows, col::KeyType)
    return table[rows][col]
end


function Base.getindex(table::DataTable, rows, cols::Array{<:KeyType,1})
    return table[rows][cols]
end


function Base.Array(table::DataTable)
    return reduce(hcat, getcolumns(table))
end


function Base.getindex(table::DataTable, rowidx::Union{Int,Colon,UnitRange{Int},Array{Int,1}}, colmap::Union{Int,Colon,UnitRange{Int},Array{Int,1}})
    table_temp = table[rowidx]
    columns = getcolumns(table_temp)[colmap]
    return reduce(hcat, columns)
end


function Base.getindex(table::DataTable, rowindex::Int, colon::Colon)
    return [ column[rowindex] for column in getcolumns(table)]
end


function Base.lastindex(table::DataTable)
    nr, nc = size(table)
    nc==0 && error("lastindex: use of 'end' in an empty DataTable")
    return nr
end


function Base.lastindex(table::DataTable, idx::Int)
    nr, nc = size(table)
    nc==0 && error("lastindex: use of 'end' in an empty DataTable")
    if idx==1
        return nr
    else
        return nc
    end
end


"""
    filter(table::DataTable, expr::Expr)

Filter rows of a DataTable object using a logical expression.
"""
function Base.filter(table::DataTable, expr::Expr)
    fields = getvars(expr)
    vars   = Dict{Symbol, Float64}()
    nr     = size(table,1)
    idxs   = falses(nr)
    for i in 1:nr
        for field in fields
            vars[field] = table[field][i]
        end
        idxs[i] = evaluate(expr; vars...)
    end
    return table[idxs]
end



"""
    split_by(table::DataTable, col_name::String)

Splits a DataTable into a vector of new DataTables based on the unique values
in a specified column.

This is an efficient implementation that makes a single pass over the data to
group row indices before creating the new tables. The resulting tables are not
guaranteed to be in any specific order.

# Arguments
- `table::DataTable`: The input table to split.
- `col_name::String`: The name of the column whose values will be used to group the rows.

# Returns
- `Vector{DataTable}`: A vector of new `DataTable` objects, where each table
  corresponds to a unique value in the `col_name` column.

# Example
```julia
header = ["ID", "Category", "Value"]
columns = [
    [1, 2, 3, 4, 5, 6],
    ["A", "B", "A", "C", "B", "A"],
    [10.1, 20.2, 10.3, 30.1, 20.4, 10.5]
]
dt = DataTable(header, columns, name="SalesData")
grouped_tables = split_by(dt, "Category")
```
"""
function split_by(table::DataTable, col_name::String)
    col_map = getcolidx(table)

    # 1. Input Validation: Check if the column exists
    if !haskey(col_map, col_name)
        throw(ArgumentError("Column '$col_name' not found in DataTable header."))
    end

    # 2. Identify the split column and its data
    col_idx   = col_map[col_name]
    columns   = getcolumns(table)
    split_col = columns[col_idx]

    # Handle empty table case
    if isempty(split_col)
        return DataTable[]
    end

    # 3. Efficiently group row indices by the unique values in the split column
    # This dictionary will map each unique value to a list of row numbers.
    # e.g., "A" => [1, 3, 6], "B" => [2, 5], ...
    groups = OrderedDict{eltype(split_col), Vector{Int}}()
    for (row_idx, value) in enumerate(split_col)
        # `get!` is a handy function: it gets the key's value, or if it doesn't
        # exist, it initializes it with the default value (the 2nd argument)
        # and then returns the value.
        indices = get!(groups, value, Int[])
        push!(indices, row_idx)
    end

    # 4. Construct the new DataTables using the grouped indices
    result_tables   = DataTable[]
    original_header = getheader(table)
    original_name   = getname(table)

    for (group_key, row_indices) in groups
        # For each column in the original table, create a new column containing
        # only the rows specified by `row_indices`. This is a fast slicing operation.
        new_cols = [col[row_indices] for col in columns]

        # Create a descriptive name for the new sub-table
        new_name = isempty(original_name) ? string(group_key) : "$(original_name)_$(group_key)"

        # Construct the new DataTable for this group and add it to the result vector
        push!(result_tables, DataTable(original_header, new_cols, name=new_name))
    end

    return result_tables
end


function Base.sort!(table::DataTable, options::NamedTuple...)
    n, m = size(table)
    idx = collect(1:n)

    for opt in options
        field = get(opt, :field, "")
        rev   = get(opt, :rev, false)

        col  = table[field][idx]
        idx_ = sortperm(col, rev=rev)
        idx  = idx[idx_]
    end

    cols = getcolumns(table)
    for i in 1:m
        cols[i] = cols[i][idx]
    end

    return table
end


"""
    compress!(table::DataTable, n::Int)

Compress a DataTable object to `n` rows. If the number of rows is less or equal to `n`, the same DataTable is returned.
"""
function compress!(table::DataTable, n::Int)
    columns = getcolumns(table)
    nr = size(table)[1]

    nr<=n && return table[:]

    factor = (nr-1)/(n-1)
    idxs   = [ round(Int, 1+factor*(i-1)) for i in 1:n ]

    for i in 1:length(columns)
        columns[i] = columns[i][idxs]
    end

    return table
end


"""
    resize(table::DataTable, n::Int=0; ratio=1.0)

    Resize a DataTable object to `n` rows. If `n` is not provided, the number of rows is calculated using `ratio` parameter.
"""
function resize(table::DataTable, n::Int=0; ratio=1.0)
    header = getheader(table)
    nr     = size(table,1)

    if n==0
        ratio > 0.0 || error("resize: ratio should be greater than zero")
        n = max(2, round(Int, nr*ratio))
    end
    n > 0 || error("resize: number of rows should be greater than zero")
    nr >= 4 || error("resize: Table object should contain at least 4 rows")

    ns = nr - 1                # number of spacings
    nb = floor(Int64, ns/3)    # number of bezier curves
    ds = 1.0 / nb              # spacing between curves

    cols = Array{Any,1}[]
    for (j,field) in enumerate(header)
        U = table[field]
        V = zeros(n)

        for (i,s) in enumerate(range(0.0, 1.0, length=n))
            # find index of Bezier and local coordinate t
            ib = floor(Int64, s/ds) + 1   # index of Bezier
            ib > nb && (ib = nb)          # fix index if s ~= 1+eps
            s0 = (ib-1) * ds              # s @ left point
            t  = (s - s0) / ds            # local t for current Bezier
            t > 1.0 && (t = 1.0)          # clean rubbish. e.g. 1.00000000002

            # collect control points
            k = 1 + (ib-1) * 3            # position of first point of bezier

            P1 = U[k  ]
            P2 = U[k+1]
            P3 = U[k+2]
            P4 = U[k+3]

            # control points
            Q1 =         P1
            Q2 = (-5.0 * P1 + 18.0 * P2 -  9.0 * P3 + 2.0 * P4) / 6.0
            Q3 = ( 2.0 * P1 -  9.0 * P2 + 18.0 * P3 - 5.0 * P4) / 6.0
            Q4 =                                            P4

            a =       Q4 - 3.0 * Q3 + 3.0 * Q2 - Q1
            b = 3.0 * Q3 - 6.0 * Q2 + 3.0 * Q1
            c = 3.0 * Q2 - 3.0 * Q1
            d =       Q1

            V[i] = a*t*t*t + b*t*t + c*t + d
        end

        push!(cols, V)
    end

    return DataTable(header, cols)
end


"""
    cut!(table::DataTable, field, value=0.0; after=false)

Cut a DataTable object at a given value of a field. If `after` is true, the DataTable is cut after the given value.
"""
function cut!(table::DataTable, field, value=0.0; after=false)
    V   = table[field] .- value
    idx = 0
    for i in 2:length(V)
        if V[i-1]*V[i] <= 0
            idx = i
            break
        end
    end

    if idx>0
        α = -V[idx-1]/(V[idx]-V[idx-1])
        header = getheader(table)
        for field in header
            W      = table[field]
            W[idx] = W[idx-1] + α*(W[idx]-W[idx-1])
        end
        rng = 1:idx
        after && (rng=1:length(V))

        columns = getcolumns(table)
        for (i,col) in enumerate(columns)
            columns[i] = col[rng]
        end
    end

    return table
end


function Base.clamp!(table::DataTable, field, lo, hi)
    clamp!(table[field], lo, hi)
    return table
end

export smooth!
function smooth!(table::DataTable, fieldx, fieldy=nothing; knots=[0.0, 1.0])
    nr = size(table,1)

    if fieldy === nothing
        fieldy = fieldx
        X = collect(range(0, 1, length=nr))
        Y = table[fieldx]
    else
        X = table[fieldx]
        Y = table[fieldy]
    end

    Idxs = split(X, X[1] .+ knots.*(X[end]-X[1]))

    for i in 1:length(Idxs)
        Xi = X[Idxs[i]]
        Yi = Y[Idxs[i]]
        M = Float64[ ones(length(Xi)) Xi Xi.^3 ]
        a, b, d = inv(M'*M)*M'*Yi
        Y[Idxs[i]] .= a .+ b*Xi .+ d*Xi.^3
    end

    table[fieldy] = Y
end


function denoise!(table::DataTable, fieldx, fieldy=nothing; noise=0.05, npatch=4)
    header = getheader(table)
    nr     = size(table,1)

    if fieldy === nothing
        X = range(0,1,length=nr)
        Y = table[fieldx]
    else
        X = table[fieldx]
        Y = table[fieldy]
    end

    # Regression along patches
    M  = ones(npatch,2)
    A  = zeros(2)
    ΔY = [ Float64[] for i in 1:nr ]

    for i in 1:nr-npatch+1
        rng = i:i+npatch-1
        Xp = X[rng]
        Yp = Y[rng]

        # Linear regression
        M[:,2] .= Xp
        A   = pinv(M)*Yp
        ΔYp = abs.(Yp .- M*A)

        for (j,k) in enumerate(rng)
            push!(ΔY[k], ΔYp[j])
        end
    end

    # Get indexes to keep
    idxs = minimum.(ΔY) .<= noise*(maximum(Y)-minimum(Y))
    # idxs[1:npatch] .= 1
    # idxs[end-npatch+1:end] .= 1

    # newtable = table[:]
    # for i in (1:nr)[idxs]
    #     j = findprev(!iszero, idxs, i-1)
    #     k = findnext(!iszero, idxs, i+1)
    #     r = (X[i]-X[j])/(X[k]-X[j])

    #     for fieldx in header
    #         V = newtable[fieldx]
    #         V[i] = V[j] + r*(V[k]-V[j])
    #     end
    # end

    # Linear interpolation of dropped points
    cols = Array{Any,1}[]
    for column in getcolumns(table)
        # V = copy(column)
        V = column
        for i in (1:nr)[.!idxs]
            j = findprev(!iszero, idxs, i-1)
            k = findnext(!iszero, idxs, i+1)
            if j===nothing
                j = k
                k = findnext(!iszero, idxs, j+1)
            end
            if k===nothing
                k = j
                j = findprev(!iszero, idxs, k-1)
            end
            r = (X[i]-X[j])/(X[k]-X[j])
            V[i] = V[j] + r*(V[k]-V[j])
        end
        push!(cols, V)
    end

    return table
    # return DataTable(header, cols)
end


# TODO: Improve column width for string items
function save(table::DataTable, filename::String; quiet=true, digits::AbstractArray=[])
    suitable_formats = (".dat", ".table", ".json", ".tex")

    _, format = splitext(filename)
    format in suitable_formats || error("DataTable: cannot save in \"$format\" format. Suitable formats $suitable_formats.")

    local f::IOStream
    try
        f  = open(filename, "w")
    catch err
        warn("DataTable: File $filename could not be opened for writing.")
        return
    end

    columns = getcolumns(table)
    header  = getheader(table)
    nr, nc  = size(table)

    if format in (".dat", ".table")
        for (i,key) in enumerate(header)
            @printf(f, "%12s", key)
            print(f, i!=nc ? "\t" : "\n")
        end

        # print values
        for i in 1:nr
            for j in 1:nc
                item = columns[j][i]
                if typeof(item)<:AbstractFloat
                    @printf(f, "%12.5e", item)
                elseif typeof(item)<:Integer
                    @printf(f, "%12d", item)
                else
                    @printf(f, "%12s", item)
                end
                print(f, j!=nc ? "\t" : "\n")
            end
        end

        quiet || printstyled("  file $filename written\n", color=:cyan)
    end

    if format==".json"
        data = OrderedDict{String,Any}()
        cols = OrderedDict{String,Any}( k=>c for (k,c) in zip(header, columns) )
        data["name"] = getname(table)
        data["data"] = cols
        JSON.print(f, data, 4)
        quiet || printstyled("  file $filename written\n", color=:cyan)
    end

    if format==".tex"
        # widths calculation
        widths = length.(header)
        types  = eltype.(columns)

        if length(digits)==0
            digits = repeat([4], nc)
        end
        @assert length(digits)==nc

        for (i,col) in enumerate(columns)
            etype = types[i]
            if etype<:AbstractFloat
                widths[i] = max(widths[i], 12)
            elseif etype<:Integer
                widths[i] = max(widths[i], 6)
            elseif etype<:AbstractString
                widths[i] = max(widths[i], maximum(length.(col)))
            else
                widths[i] = max(widths[i], maximum(length.(string.(col))))
            end
        end

        # printing header
        level = 1
        indent = "    "
        println(f, indent^level, raw"\begin{tabular}{", "c"^nc, "}" )
        level = 2
        println(f, indent^level, raw"\toprule")
        print(f, indent^level)
        for (i,key) in enumerate(header)
            etype = types[i]
            width = widths[i]
            if etype<:Real
                print(f, lpad(key, width))
            else
                print(f, rpad(key, width))
            end
            i<nc && print(f, " & ")
        end
        println(f, raw" \\\\")

        # printing body
        println(f, indent^level, raw"\hline")
        for i in 1:nr
            print(f, indent^level)
            for j in 1:nc
                etype = types[j]
                item = columns[j][i]
                width = widths[j]
                if etype<:AbstractFloat
                    #item = @sprintf("%12.3f", item)
                    # Use Printf.format for dynamic format strings
                    dig = digits[j]
                    if isnan(item)
                        item = "-"
                    else
                        fmt = Printf.Format("%$(width).$(dig)f")
                        item = Printf.format(fmt, item)
                    end
                    print(f, lpad(string(item), width))
                elseif etype<:Integer
                    item = @sprintf("%6d", item)
                    print(f, lpad(item, width))
                elseif etype<:AbstractString
                    print(f,rpad(item, width))
                else
                    str = string(item)
                    print(f, rpad(item, width))
                end
                j<nc && print(f, " & ")
            end
            println(f, raw" \\\\")
        end
        println(f, indent^level, raw"\bottomrule")

        # printing ending
        level = 1
        println(f, indent^level, raw"\end{tabular}")
    end

    close(f)
    return nothing
end


function DataTable(filename::String, delim::Char='\t')
    _, format = splitext(filename)
    formats = (".dat", ".table", ".csv")
    format in formats || error("DataTable: cannot read \"$format\". Suitable formats are $formats")

    if format in formats
        matrix, headstr = readdlm(filename, delim, header=true, use_mmap=false)
        header = vec(strip.(headstr))
        table = DataTable(header, matrix)
        return table
    end
end


# TODO: Improve display. Include columns datatypes
function Base.show(io::IO, table::DataTable)
    columns = getcolumns(table)

    println(io)
    nr, nc = size(table)
    if nr==0
        print(io, "DataTable()")
        return
    end

    print(io, "DataTable ")
    name = getname(table)
    if name != ""
        println(io, repr(name))
    else
        println(io)
    end

    header = getheader(table)
    types = eltype.(columns)

    hwidths = length.(header)
    widths  = zeros(Int, length(header))
    useformat   = falses(length(header))
    shortformat = falses(length(header))

    if nr>0
        for (i,col) in enumerate(columns)
            etype = types[i]
            widths[i] = maximum(length.(string.(col)))
            if etype<:AbstractFloat
                if widths[i] >= 11
                    widths[i] = 11
                    useformat[i] = true
                end
            end
            widths[i] = max(widths[i], hwidths[i])
        end
    end

    total = sum(widths) + (length(header)+1)*3
    if total>displaysize(stdout)[2]
        for (i,col) in enumerate(columns)
            etype = types[i]
            if etype<:AbstractFloat && useformat[i] && hwidths[i]<=8
                shortformat[i] = true
                widths[i] = max(8, hwidths[i])
            end
        end
    end

    # print header
    print(io, " │ ")
    for (i,key) in enumerate(header)
        etype = types[i]
        width = widths[i]
        if etype<:Real
            print(io, lpad(key, width))
        else
            print(io, rpad(key, width))
        end
        print(io, " │ ")
    end
    println(io)

    visible_rows = 30
    half_vrows = div(visible_rows,2)

    # print values
    for i in 1:nr
        if i>half_vrows && nr-i>=half_vrows
            i==half_vrows+1 && println(io, " ⋮")
            continue
        end

        print(io, " │ ")
        for j in 1:nc
            etype = types[j]
            item = columns[j][i]
            if etype<:AbstractFloat
                if useformat[j]
                    if shortformat[j]
                        item = @sprintf("%8.1e", item)
                    else
                        item = @sprintf("%11.4e", item)
                    end
                end
                print(io, lpad(item, widths[j]))
            else
                print(io, rpad(item, widths[j]))
            end
            print(io, " │ ")
        end
        i<nr && println(io)
    end

end


randtable() = DataTable(["x","y","z"], Any[0:10 rand().*(sin.(0:10).+(0:10)) rand().*(cos.(0:10).+(0:10)) ])
