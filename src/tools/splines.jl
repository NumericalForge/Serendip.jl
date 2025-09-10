# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export LinearSpline, CubicSpline

abstract type AbstractSpline end

function check_Spline(X::Vector{<:Real}, Y::Vector{<:Real}; left::Symbol=:flat, right::Symbol=:flat)
    length(X) == length(Y) || error("Spline: X and Y must have the same length.")
    length(X) >= 2 || error("Spline: At least two points are required.")
    issorted(X) || error("Spline: X must be sorted in ascending order.")
    left in (:flat, :linear, :error) || error("Spline: left extrapolation must be :flat, :linear, or :error.")
    right in (:flat, :linear, :error) || error("Spline: right extrapolation must be :flat, :linear, or :error.")
end


# ❱❱❱ Linear Spline Implementation

struct LinearSpline<:AbstractSpline
    X::Vector{Float64}  # x values
    Y::Vector{Float64}  # y values
    left::Symbol    # left extrapolation method (:flat, :linear, :error)
    right::Symbol   # right extrapolation method (:flat, :linear, :error)
    
    function LinearSpline(X::Vector{<:Real}, Y::Vector{<:Real}; left::Symbol=:flat, right::Symbol=:flat)
        check_Spline(X, Y; left=left, right=right)
        return new(X, Y, left, right)
    end
    
    function LinearSpline(XY::Matrix{<:Real}; left::Symbol=:flat, right::Symbol=:flat)
        @check size(XY, 2) == 2 "LinearSpline: Input matrix must have two columns."
        X = XY[:, 1]
        Y = XY[:, 2]
        check_Spline(X, Y; left=left, right=right)
        return LinearSpline(X, Y, left=left, right=right)
    end
end


function (spline::LinearSpline)(x::Real)
    X = spline.X
    Y = spline.Y
    n = length(X)

    # ❱❱❱ Handle Extrapolation
    if x < X[1]
        if spline.left == :error
            throw(DomainError(x, "Value is below the interpolation range."))
        elseif spline.left == :flat
            return Y[1]
        elseif spline.left == :linear
            # Extrapolate using the slope of the first segment
            slope = (Y[2] - Y[1]) / (X[2] - X[1])
            return Y[1] + slope * (x - X[1])
        end
    elseif x > X[n]
        if spline.right == :error
            throw(DomainError(x, "Value is above the interpolation range."))
        elseif spline.right == :flat
            return Y[n]
        elseif spline.right == :line
            # Extrapolate using the slope of the last segment
            slope = (Y[n] - Y[n-1]) / (X[n] - X[n-1])
            return Y[n] + slope * (x - X[n])
        end
    end

    # ❱❱❱ Find the Interval
    # `searchsortedlast` is a fast binary search to find the index of the segment
    # such that X[i] <= x <= X[i+1]
    i = searchsortedlast(X, x)
    i==n && return Y[n] # handle the case when x == X[n]
    
    # ❱❱❱ Linear Interpolation
    x1, y1 = X[i], Y[i]
    x2, y2 = X[i+1], Y[i+1]

    # Handle vertical segment to avoid division by zero
    if x1 == x2
        return y1
    end

    # Apply the linear interpolation formula
    return y1 + (x - x1) * (y2 - y1) / (x2 - x1)
end


function derive(spline::LinearSpline, x::Real)
    X = spline.X
    Y = spline.Y
    n = length(X)

    # ❱❱❱ Handle Extrapolation
    if x < X[1]
        if spline.left == :error
            throw(DomainError(x, "Value is below the interpolation range."))
        elseif spline.left == :flat
            return 0.0
        elseif spline.left == :linear
            return (Y[2] - Y[1]) / (X[2] - X[1])
        end
    elseif x > X[n]
        if spline.right == :error
            throw(DomainError(x, "Value is above the interpolation range."))
        elseif spline.right == :flat
            return 0.0
        elseif spline.right == :linear
            return (Y[n] - Y[n-1]) / (X[n] - X[n-1])
        end
    end

    # ❱❱❱ Find the Interval
    i = searchsortedlast(X, x)

    # If x is exactly the last point, use the last segment's slope
    if i == n
        i = n - 1
    end

    # ❱❱❱ Calculate Slope of the Segment
    x1, y1 = X[i], Y[i]
    x2, y2 = X[i+1], Y[i+1]

    if x1 == x2
        return Inf # Derivative of a vertical line is infinite
    end

    return (y2 - y1) / (x2 - x1)
end


# ❱❱❱ Cubic Spline Implementation

struct CubicSpline<:AbstractSpline
    X::Vector{Float64}
    C::Vector{Float64}
    left::Symbol
    right::Symbol
    
    function CubicSpline(X::Vector{<:Real}, Y::Vector{<:Real}; left::Symbol=:flat, right::Symbol=:flat)
        check_Spline(X, Y; left=left, right=right)

        n = length(X)
        h = diff(X)

        @show h

        @check !any(h .== 0) "Spline: All X values must be distinct."
        
        M = zeros(n) # Initialize second derivatives vector

        # The tridiagonal system is only needed for 3 or more points
        if n > 2
            # Set up and solve the tridiagonal system for the second derivatives (M)
            dl = h[2:n-2]
            d = 2 * (h[1:n-2] + h[2:n-1])
            du = h[2:n-2] # superdiagonal is the same as subdiagonal in this system
            A = Tridiagonal(dl, d, du)
            
            # Right-hand side vector, calculated robustly
            B = 6 * diff(diff(Y) ./ h)
            
            # Solve for interior second derivatives
            M_interior = A \ B
            M[2:n-1] = M_interior
            # For a natural spline, M[1] and M[n] remain 0
        end
        # If n=2, M is just [0.0, 0.0], which is correct for a line

        # Compute and store coefficients for each segment
        coeffs = zeros(4 * (n - 1))
        for i in 1:(n-1)
            a = Y[i]
            b = (Y[i+1] - Y[i]) / h[i] - h[i] * (2 * M[i] + M[i+1]) / 6
            c = M[i] / 2
            d = (M[i+1] - M[i]) / (6 * h[i])
            
            offset = 4 * (i - 1)
            coeffs[offset+1] = a
            coeffs[offset+2] = b
            coeffs[offset+3] = c
            coeffs[offset+4] = d
        end
        
        return new(X, coeffs, left, right)
    end
    
    function CubicSpline(XY::Matrix{<:Real}; left::Symbol=:flat, right::Symbol=:flat)
        @check size(XY, 2) == 2 "CubicSpline: Input matrix must have two columns."

        X = XY[:, 1]
        Y = XY[:, 2]
        return CubicSpline(X, Y; left=left, right=right)
    end
end


function (spline::CubicSpline)(x::Real)
    X = spline.X
    n = length(X)
    
    # ❱❱❱ Handle Extrapolation
    if x < X[1]
        if spline.left == :error
            throw(DomainError(x, "Value is below the interpolation range."))
        elseif spline.left == :flat
            return spline.C[1] # Y[1] is the first 'a' coefficient
        elseif spline.left == :linear
            y1 = spline.C[1]
            slope = derive(spline, X[1])
            return y1 + slope * (x - X[1])
        end
    elseif x > X[n]
        if spline.right == :error
            throw(DomainError(x, "Value is above the interpolation range."))
        elseif spline.right == :flat
            # The last y-value is the 'a' coefficient of the last segment plus its value at h
            i = n - 1
            offset = 4 * (i - 1)
            a, b, c, d = spline.C[offset+1:offset+4]
            h = X[n] - X[n-1]
            return a + b*h + c*h^2 + d*h^3 # Y[n]
        elseif spline.right == :linear
            y_n = spline(X[n]) # Get Y[n] value
            slope = derive(spline, X[n])
            return y_n + slope * (x - X[n])
        end
    end

    # ❱❱❱ Find the Interval and Interpolate
    i = searchsortedlast(X, x)
    if i == n; i = n - 1; end

    offset = 4 * (i - 1)
    a, b, c, d = spline.C[offset+1 : offset+4]
    
    dx = x - X[i]
    
    # Evaluate polynomial using Horner's method for efficiency
    return a + dx * (b + dx * (c + dx * d))
end


function derive(spline::CubicSpline, x::Real)
    X = spline.X
    n = length(X)

    # ❱❱❱ Handle Extrapolation
    if x < X[1]
        if spline.left == :error
            throw(DomainError(x, "Value is below the interpolation range."))
        elseif spline.left == :flat
            return 0.0
        elseif spline.left == :linear
            # The derivative is constant for linear extrapolation
            return derive(spline, X[1])
        end
    elseif x > X[n]
        if spline.right == :error
            throw(DomainError(x, "Value is above the interpolation range."))
        elseif spline.right == :flat
            return 0.0
        elseif spline.right == :linear
            return derive(spline, X[n])
        end
    end

    # ❱❱❱ Find the Interval and Differentiate
    i = searchsortedlast(X, x)
    if i == n; i = n - 1; end

    offset = 4 * (i - 1)
    _, b, c, d = spline.C[offset+1 : offset+4]

    dx = x - X[i]
    
    # Derivative of a+b*dx+c*dx^2+d*dx^3 is b + 2c*dx + 3d*dx^2
    return b + dx * (2c + dx * 3d)
end