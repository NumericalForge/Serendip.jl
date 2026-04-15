# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


function domain_plot_view_basis(up::Symbol)
    ex = Vec3(1.0, 0.0, 0.0)
    ey = Vec3(0.0, 1.0, 0.0)
    ez = Vec3(0.0, 0.0, 1.0)

    if up == :z
        return ex, ey, ez
    elseif up == :x
        return ey, ez, ex
    elseif up == :y
        return ez, ex, ey
    end

    error("domain_plot_view_basis: unsupported up direction $up")
end


function view_to_world(X; up::Symbol=:z)
    X = Vec3(X)
    forward, right, upvec = domain_plot_view_basis(up)
    return X[1]*forward + X[2]*right + X[3]*upvec
end


function world_to_view(X; up::Symbol=:z)
    X = Vec3(X)
    forward, right, upvec = domain_plot_view_basis(up)
    return Vec3(dot(X, forward), dot(X, right), dot(X, upvec))
end


function rotate_view_point(X, azimuth, elevation)
    X = Vec3(X)

    θ = -azimuth*pi/180
    R = Quaternion(cos(θ/2), 0, 0, sin(θ/2))
    X = Vec3((R*X*conj(R))[2:4])

    θ = elevation*pi/180
    R = Quaternion(cos(θ/2), 0, sin(θ/2), 0)
    return Vec3((R*X*conj(R))[2:4])
end


function project_view_point(X, azimuth, elevation, distance; up::Symbol=:z, center=Vec3(0.0, 0.0, 0.0), reflength::Real=1.0)
    X = world_to_view(Vec3(X) - Vec3(center), up=up)
    X = rotate_view_point(X, azimuth, elevation)

    distance == 0 && (distance = reflength*3)
    distance = max(distance, reflength)
    focal_length = 0.1*distance

    x = X[1]
    depth = distance - x
    y = X[2]*focal_length/depth
    z = X[3]*focal_length/depth

    return Vec3(y, z, depth)
end
