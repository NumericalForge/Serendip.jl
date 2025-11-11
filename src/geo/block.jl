# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

function box_coords(C1::AbstractArray{<:Real}, C2::AbstractArray{<:Real}, ndim::Int)
    x1 = C1[1]
    y1 = C1[2]
    lx = C2[1] - C1[1]
    ly = C2[2] - C1[2]

    if ndim==2
        return [
                 x1     y1     0.0
                 x1+lx  y1     0.0
                 x1+lx  y1+ly  0.0
                 x1     y1+ly  0.0
               ]
    else
        z1 = C1[3]
        lz = C2[3] - C1[3]
        return [
                 x1     y1     z1
                 x1+lx  y1     z1
                 x1+lx  y1+ly  z1
                 x1     y1+ly  z1
                 x1     y1     z1+lz
                 x1+lx  y1     z1+lz
                 x1+lx  y1+ly  z1+lz
                 x1     y1+ly  z1+lz
                ]
    end
end


mutable struct Block
    ndim::Int
    points::Vector{Point}
    blockshape::CellShape
    shape::CellShape
    nx::Int64
    ny::Int64
    nz::Int64
    rx::Float64
    ry::Float64
    rz::Float64
    tag::String
    id::Int64



    function Block(
        points::Vector{Point};
        nx::Int  = 0,
        ny::Int  = 0,
        nz::Int  = 0,
        n ::Int  = 0,
        rx::Real = 1.0,
        ry::Real = 1.0,
        rz::Real = 1.0,
        r ::Real = 0.0,
        shape = nothing,
        tag       = "",
        )

        (nx>0 || n>0 ) || error("Block: nx or n must be greater than zero.")
        r>0 && (rx=r)

        shapes1d = (LIN2, LIN3, LIN4)
        shapes2d = (TRI3, TRI6, QUAD4, QUAD8, QUAD9, QUAD12)
        shapes3d = (TET4, TET10, HEX8, HEX20, HEX27, PYR5)

        # Get ndim
        sumy = sum(abs, [ p.coord.y for p in points ])
        sumz = sum(abs, [ p.coord.z for p in points ])

        ndim = 3
        sumz==0 && (ndim=2)
        sumy+sumz==0 && (ndim=1)

        if n>0
            nx = n
            ndim = 1
        end

        # Check for surface or chord
        # surface = ndim==3 && nz==0
        # chord   = ndim>1 && ny==0 && nz==0

        nz==0 && ndim==3 && (nz=1)
        ny==0 && ndim>=2 && (ny=1)
        shape in shapes3d && (ndim==3 || error("Block: 3d points and nx, ny and nz are required for cell blockshape $(shape.name)"))

        # prisma points if given just two points
        # if ndim in (2,3) && length(points)==2 && !surface && !chord
        #     coords = box_coords(points[1].coord, points[2].coord, ndim)
        #     points = Point[]
        #     for i in 1:size(coords, 1)
        #         push!(points, Point(coords[i, :]))
        #     end
        # end

        npoints = length(points)

        if ndim==1
            npoints in (2, 3) || error("Block: invalid number of points ($npoints) for dimension $ndim or chord.")
            shape===nothing && (shape=LIN2)
            shape in shapes1d || error("Block: invalid cell type $(shape.name) for dimension $ndim or chord.")
            blockshape = npoints==2 ? LIN2 : LIN3
        elseif ndim==2
            npoints in (4, 8) || error("Block: invalid number of points ($npoints) for dimension $ndim or surface.")
            shape===nothing && (shape=QUAD4)
            shape in shapes2d || error("Block: invalid cell type $(shape.name) for dimension $ndim or surface.")
            blockshape = npoints==4 ? QUAD4 : QUAD8
        else
            npoints in (8, 20) || error("Block: invalid number of points ($npoints) for dimension $ndim.")
            shape===nothing && (shape=HEX8)
            shape in shapes3d || error("Block: invalid cell type $(shape.name) for dimension $ndim.")
            blockshape = npoints==8 ? HEX8 : HEX20
        end

        for i in 1:length(points)
            points[i].id = i
        end

        return new(ndim, points, blockshape, shape, nx, ny, nz, rx, ry, rz, tag)
    end


    # function Block(X1::Vector{<:Real}, X2::Vector{<:Real}; args...)
    #     nz = get(args, :nz, 0)
    #     ny = get(args, :ny, 0)

    #     length(X1) == length(X2) || error("Block: X1 and X2 must have the same length")

    #     ndim = max(length(X1))
    #     ndim>3 && error("Block: invalid dimension $ndim for block")

    #     # Check for surface or chord
    #     surface = ndim==3 && nz==0
    #     chord   = ndim>1 && ny==0 && nz==0

    #     if ndim in (2,3) && !surface && !chord
    #         coords = box_coords(X1, X2, ndim)
    #         ncoord = size(coords,1)
    #     else
    #         coords = [ X1'; X2' ]
    #         ncoord = 2
    #     end

    #     points = [ Point(coords[i,:]) for i in 1:ncoord ]

    #     return Block(points; args...)
    # end


    function Block(X::Vector{<:Real}, dx::Float64=0.0, dy::Float64=0.0, dz::Float64=0.0; args...)

        # nz = get(args, :nz, 0)
        # ny = get(args, :ny, 0)

        # dx = get(args, :dx, 0.0)
        # dy = get(args, :dy, 0.0)
        # dz = get(args, :dz, 0.0)
        n = get(args, :n, 0)

        # ndim = sign(dx)^2 + sign(dy)^2 + sign(dz)^2
        ndim = 1
        dz!=0.0 && (ndim=3)
        dz==0.0 && dy!=0.0 && (ndim=2)
        dz==0.0 && dy==0.0 && (ndim=1)

        if dz==0.0
            if dy==0.0
                ndim = 1
            else
                ndim = 2
            end
        else
            ndim = 3
        end

        if n>0
            ndim = 1
            # dl>0 || error("Block: dl must be greater than zero when n is specified.")
        end

        if length(X)!=3
            X = vcat(X, zeros(3 - length(X)))
        end
        X2 = X + [dx, dy, dz]
        if ndim in (2,3)
            coords = box_coords(X, X2, ndim)
            ncoord = size(coords,1)
        else
            coords = [ X'; X2' ]
            ncoord = 2
        end

        points = [ Point(coords[i,:]) for i in 1:ncoord ]

        return Block(points; args...)
    end


    function Block(coords::Matrix{<:Real}; args...)
        # nz = get(args, :nz, 0)
        # ny = get(args, :ny, 0)

        ncoord, ncol = size(coords)
        ncol<=3 || error("Block: invalid coordinate matrix")
        ncoord in (2, 3, 4, 8, 20) || error("Block: invalid number of points ($ncoord) in coordinate matrix")

        # Get ndim
        # sumy = ncol>=2 ? sum(abs, coords[:,2]) : 0.0
        # sumz = ncol==3 ? sum(abs, coords[:,3]) : 0.0

        # ndim = 3
        # sumz==0 && (ndim=2)
        # sumy+sumz==0 && (ndim=1)


        # # Check for surface or chord
        # surface = ndim==3 && nz==0
        # chord   = ndim>1 && ny==0 && nz==0

        # if ndim in (2,3) && ncoord==2 && !surface && !chord
        #     coords = box_coords(coords[1,:], coords[2,:], ndim)
        #     ncoord = size(coords,1)
        # end
        points = [ Point(coords[i,:]) for i in 1:ncoord ]

        return Block(points; args...)
    end
end


function Base.copy(block::Block)

    return Block(copy(get_coords(block.points)), nx=block.nx, ny=block.ny, nz=block.nz, shape=block.shape, tag=block.tag)
    # return Block(copy(block.points), nx=block.nx, ny=block.ny, nz=block.nz, shape=block.shape, tag=block.tag)
end



function Base.copy(blocks::Vector{Block})
    return [ copy(bl) for bl in blocks ]
end


function BlockGrid(
    X::Array{<:Real},            # list of x coordinates
    Y::Array{<:Real},            # list of y coordinates
    Z::Array{<:Real}=Float64[];  # list of z coordinates
    nx=[],  # list of divisions in the x direction
    ny=[],  # list of divisions in the y direction
    nz=[],  # list of divisions in the z direction
    rx=[],  # list of divisions ratios in the x direction
    ry=[],  # list of divisions ratios in the x direction
    rz=[],  # list of divisions ratios in the x direction
    shape=QUAD4, # element blockshape
    tag="",          # elements tag
    id=-1
    )

    length(rx)==0 && (rx = ones(length(nx)))
    length(ry)==0 && (ry = ones(length(ny)))
    length(rz)==0 && (rz = ones(length(nz)))
    blocks = Block[]
    for i in 1:length(nx)
        for j in 1:length(ny)
            coords = [
                X[i] Y[j]
                X[i+1] Y[j+1]
            ]
            bl = Block(coords, shape=shape, nx=nx[i], ny=ny[j], rx=rx[i], ry=ry[j])
            push!(blocks, bl)
        end
    end
    return blocks
end

