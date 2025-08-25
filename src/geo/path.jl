export Path, addinset!

mutable struct PathCmd
    key::Symbol
    idxs::Array{Int}
    C::Vector{Float64} # center for arcs
    N::Vector{Float64} # normal for arcs

    function PathCmd(key::Symbol, idxs::Vector{Int})
        return new(key, idxs, [], [])
    end
end

function Base.copy(cmd::PathCmd)
    new_cmd = PathCmd(cmd.key, copy(cmd.idxs))
    new_cmd.C = cmd.C
    new_cmd.N = cmd.N
    return new_cmd
end


mutable struct Path
    points::Array{Point}
    cmds::Array{PathCmd}
    closed::Bool
    len::Array{Float64,1} # normalized cumulative length

    function Path(points::Vector{Point}, cmds::Vector{PathCmd}; closed::Bool=false)
        # check commands
        for cmd in cmds
            if cmd.key in (:A, :Ac)
                length(cmd.idxs) != 3 && error("PathCmd: Arc command requires 3 points")

                X1, X2, X3 = [ p.coord for p in points[cmd.idxs] ]
                if cmd.key == :Ac # point-center-point arc
                    N = normalize(cross(X1-X2, X3-X2))
                    cmd.C = X2
                    cmd.N = N
                else # 3 points arc
                    V1 = X2 - X1
                    V2 = X3 - X2
                    A = 0.5*norm(cross(V1, V2)) # area of the triangle
                    A==0 && error("Invalid arc path: points are collinear")
                    # escalars
                    b = 8*A^2
                    k1 = b*dot(V2,V2)*(dot(V1,V1) - dot(V1,V2))
                    k2 = b*dot(V1,V1)*(dot(V2,V2) - dot(V1,V2))
                    C = X1 + k1*V1 + k2*V2 # center
                    N = normalize(cross(X1-C, X3-C))
                    cmd.C = C
                    cmd.N = N
                end
            elseif cmd.key == :B
                length(cmd.idxs) in (3, 4) || error("PathCmd: Bezier command requires 3 or 4 points")
            end
        end

        this = new(points, cmds, closed, [])
        l = 0
        len = Float64[]
        for cmd in cmds
            l += length(this, cmd)
            push!(this.len, l)
        end
        this.len = this.len ./ l # normalize

        return this
    end
end


function Path(pts::Vector{<:Real}...; closed::Bool=false)
    points = [ Point(pts[i]) for i in 1:length(pts) ]
    cmds = [ PathCmd(:M, points[1]) ]
    for i in 2:length(points)
        push!(cmds, PathCmd(:L, points[i]))
    end
    # create a new Path
    return Path(points, cmds; closed=closed)
end


function Base.copy(path::Path)
    cmds = copy.(path.cmds) # also copies the points
    points = copy.(path.points)
    # points = Point[]
    # for cmd in cmds
    #     if cmd.key==:M
    #         push!(points, cmd.points[1])
    #     else
    #         cmd.points[1] = points[end]
    #         append!(points, cmd.points[2:end])
    #     end
    # end

    return Path(points, cmds, closed=path.closed)

end


function Path(tokens::Union{Symbol,Number}...; closed=false)
    # any( t->isa(t, Point), tokens ) || return path_from_numbers(tokens...; closed=closed)

    n  = length(tokens)
    idx = 1
    cmds = PathCmd[]
    points = Point[]

    # local startpoint, endpoint
    while idx<=n
        token = tokens[idx]
        if token==:M
            p1 = tokens[idx+1]
            cmd = PathCmd(token, [ p1 ])
            push!(points, p1)
            push!(cmds, cmd)
            idx+=2
        elseif token==:L
            p1 = points[end]
            p2 = tokens[idx+1]
            cmd = PathCmd(token, [ p1, p2 ])
            push!(points, p2)
            push!(cmds, cmd)
            idx+=2
        elseif token==:A
            p1 = points[end]
            p2 = tokens[idx+1]
            p3 = tokens[idx+2]
            cmd = PathCmd(token, [p1, p2, p3])
            push!(points, p2)
            push!(points, p3)
            push!(cmds, cmd)
            idx+=3
        elseif token==:C
            p1 = points[end]
            p2 = tokens[idx+1]
            p3 = tokens[idx+2]
            p4 = tokens[idx+3]
            cmd = PathCmd(token, [p1, p2, p3, p4])
            push!(points, p2)
            push!(points, p3)
            push!(points, p4)
            push!(cmds, cmd)
            idx+=4
        else
            error("Invalid token $token")
        end

    end

    if closed && points[1]!==points[end]
        cmd = LineCmd(endpoint, startpoint)
        push!(cmds, cmd)
    end

    if !closed && points[1]===points[end]
        closed = true
    end

    # get unique points in case of repeated points
    points = unique(points)

    # normalized length
    L = cumsum( length(pc) for pc in cmds )
    L = L./L[end] # normalize

    return Path(points, cmds, closed, L)
end


function Base.length(path::Path, cmd::PathCmd)
    points = path.points

    if cmd.key==:M
        return 0.0
    elseif cmd.key==:L
        X1, X2 = [ p.coord for p in points[cmd.idxs] ]
        return norm(X2 - X1)
    elseif cmd.key==:Ac # point-center-point arc
        X1, X2, X3 = [ p.coord for p in points[cmd.idxs] ]
        r = norm(X2-X1)
        θ = acos(clamp(dot(X1-X2, X3-X2)/r^2, -1, 1))
        len = r*θ
        return len
    elseif cmd.key==:A # 3 points arc
        X1, X2, X3 = [ p.coord for p in points[cmd.idxs] ]
        C = cmd.C
        r = norm(C - X1)
        θ = acos(clamp(dot(X1-C, X3-C)/r^2, -1, 1))
        len = r*θ
        return len
    elseif cmd.key==:B
        ipoints = [0.5 - √(3/5)/2, 0.5, 0.5 + √(3/5)/2]
        weights = [5/18, 8/18, 5/18]
        if length(cmd.idxs) == 3 # quadratic Bezier
            X1, X2, X3 = [ p.coord for p in points[cmd.idxs] ]
            deriv3(t) = 2*(1 - t)*(X2 - X1) + 2*t*(X3 - X2)
            # Approximate the arc length
            return sum([weights[i] * norm(deriv3(ipoints[i])) for i in 1:3])
        elseif length(cmd.idxs) == 3 # quadratic Bezier
            X1, X2, X3, X4 = [ p.coord for p in points[cmd.idxs] ]
            deriv4(t) = 3*(1 - t)^2*(X2 - X1) + 6*(1 - t)*t*(X3 - X2) + 3*t^2*(X4 - X3)
            return sum([weights[i]*norm(deriv4(ipoints[i])) for i in 1:3])
        else
            error("Path: Unsupported Bezier with $(length(cmd.idxs)) points")
        end
    else
        error("Invalid command key $(cmd.key)")
    end

end


function evaluate(path::Path, cmd::PathCmd, t::Float64)
    points = path.points

    if cmd.key==:M
        return points[cmd.idxs[1]].coord
    elseif cmd.key==:L
        X1, X2 = [ p.coord for p in points[cmd.idxs] ]
        return X1 + t*(X2 - X1)
    elseif cmd.key==:Ac  # Arc path
        X1, X2, X3 = [ p.coord for p in points[cmd.idxs] ]
        N = cmd.N
        norm(N) == 0.0 && error("Invalid arc path: points are collinear")
        r = norm(X1-X2)
        θ = t*acos(clamp(dot(X1-X2, X3-X2)/r^2, -1, 1))
        R = Quaternion(cos(θ/2), N[1]*sin(θ/2), N[2]*sin(θ/2), N[3]*sin(θ/2))
        P = X2 + R*(X1 - X2)*conj(R)
        return P
    elseif cmd.key==:A # 3 points arc
        X1, X2, X3 = [ p.coord for p in points[cmd.idxs] ]
        C, N = cmd.C, cmd.N
        θ = t*acos(clamp(dot(X1-C, X3-C)/r^2, -1, 1))
        R = Quaternion(cos(θ/2), N[1]*sin(θ/2), N[2]*sin(θ/2), N[3]*sin(θ/2))
        P = C + R*(X1 - C)*conj(R)
        return P
    elseif cmd.key==:B
        if length(cmd.idxs) == 3 # quadratic Bezier
            X1, X2, X3 = [ p.coord for p in points[cmd.idxs] ]
            return (1-t)^2*X1 + 2*(1-t)*t*X2 + t^2*X3
        elseif length(cmd.idxs) == 4 # cubic Bezier
            X1, X2, X3, X4 = [ p.coord for p in points[cmd.idxs] ]
            return (1-t)^3*X1 + 3*(1-t)^2*t*X2 + 3*(1-t)*t^2*X3 + t^3*X4
        end
    else
        error("Invalid command key $(cmd.key)")
    end
end


function evaluate(path::Path, t::Float64)
    i = searchsortedfirst(path.len, t)
    if i==1
        return path.cmds[1](0.0)
    elseif i>length(path.len)
        return path.cmds[end](1.0)
    else
        tt = (t-path.len[i-1])/(path.len[i]-path.len[i-1])
        return evaluate(path, path.cmds[i], tt)
    end
end

