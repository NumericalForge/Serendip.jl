# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


function setup_tensile_strength(ft::Real, ft_rule::Union{Symbol,AbstractSpline}, GF::Real, wc::Real)
    # ft_rule: :constant, :hordijk, :linear, :bilinear, :soft, or custom function (Spline)
    # GF: Fracture energy
    # wc: critical crack opening

    (ft_rule in (:linear, :bilinear, :hordijk, :soft) || ft_rule isa AbstractSpline) ||
        return wc, ft_rule, nothing, failure("ft_rule must be :linear, :bilinear, :hordijk, :soft, or a custom function (Spline). Got $(repr(ft_rule)).")

    ft_fun = nothing
    if ft_rule isa AbstractSpline
        ft_fun = ft_rule
        ft_rule = :custom

        x0, y0 = ft_fun.X[1], ft_fun.Y[1]
        ft_fun.X[1]!=0.0 && return wc, ft_rule, ft_fun, failure("ft_fun must start at (0, tensile_strength=$ft). Got ($x0, $y0)")
        abs(ft_fun.Y[1]-ft)>1e-4 && return wc, ft_rule, nothing, failure("ft_fun must start at (0, tensile_strength=$ft). Got ($x0, $y0)")

        if ft_fun.Y[end] == 0.0
            wc = ft_fun.X[end]
        else
            wc = Inf
        end
    else
        if isnan(wc)
            isnan(GF) && failure("wc or GF must be defined when using a predefined ft_rule model")
            @check GF>0 "GF must be positive. Got GF=$GF"

            if ft_rule == :constant
                wc = Inf
            elseif ft_rule==:hordijk
                wc = GF/(0.1947*ft)
            elseif ft_rule == :linear
                wc = 2*GF/ft
            elseif ft_rule == :bilinear
                wc = 5*GF/ft
            elseif ft_rule==:soft
                wc = GF/(0.1947*ft)
            end

            wc = round(wc, sigdigits=5)
        end
    end

    wc > 0.0 || return wc, ft_rule, nothing, failure("wc must be positive. Got wc=$(repr(wc))")

    return wc, ft_rule, ft_fun, success()
end


function calc_tensile_strength(ft::Real, ft_rule::Symbol, ft_fun::Union{AbstractSpline,Nothing}, wc::Real, w::Real)
    # ft: tensile strength at zero crack opening
    # ft_rule: :constant, :hordijk, :linear, :bilinear, :soft, or custom function (Spline)
    # ft_fun: custom function (Spline)
    # w: crack opening (plastic displacement jump)

    if ft_rule == :constant
        return ft
    elseif ft_rule == :hordijk
        if w < wc
            e = exp(1.0)
            z = (1 + 27*(w/wc)^3)*e^(-6.93*w/wc) - 28*(w/wc)*e^(-6.93)           
        else
            z = 0.0
        end
        return z*ft 
    elseif ft_rule == :linear
        if w < wc
            a = ft 
            b = ft/wc
        else
            a = 0.0
            b = 0.0
        end
        return a - b*w
    elseif ft_rule == :bilinear
        σs = 0.25*ft 
        if w < ws
            a  = ft  
            b  = (ft  - σs)/ws
        elseif w < wc
            a  = wc*σs/(wc-ws)
            b  = σs/(wc-ws)
        else
            a = 0.0
            b = 0.0
        end
        return a - b*w
    elseif ft_rule == :soft
        m = 0.55
        a = 1.30837
        if w == 0.0
            z = 1.0
        elseif 0.0 < w < wc
            x = w/wc
            z = 1.0 - a^(1.0 - 1.0/x^m)
        else
            z = 0.0
        end
        return z*ft
    else
        @assert ft_fun isa AbstractSpline
        return ft_fun(w)
    end

    return ft
end


function calc_tensile_strength_derivative(ft::Real, ft_rule::Symbol, ft_fun::Union{AbstractSpline,Nothing}, wc::Real, w::Real)
    # ∂σmax/∂up = dft_dw
    if ft_rule == :linear
        if w < wc
            b = ft /wc
        else
            b = 0.0
        end
        dft_dw = -b
    elseif ft_rule == :bilinear
        σs = 0.25*ft 
        if w < ws
            b  = (ft  - σs)/ws
        elseif w < wc
            b  = σs/(wc-ws)
        else
            b = 0.0
        end
        dft_dw = -b
    elseif ft_rule == :hordijk
        if w < wc
            e = exp(1.0)
            dz = ((81*w^2*e^(-6.93*w/wc)/wc^3) - (6.93*(1 + 27*w^3/wc^3)*e^(-6.93*w/wc)/wc) - 0.02738402432/wc)
        else
            dz = 0.0
        end
        dft_dw = dz*ft 
    elseif ft_rule == :soft
        m = 0.55
        a = 1.30837

        if w == 0.0
            dz = 0.0
        elseif w < wc
            x = w/wc
            dz =  -m*log(a)*a^(1-x^-m)*x^(-m-1)/wc
        else
            dz = 0.0
        end
        dft_dw = dz*ft 
    else
        dft_dw = derive(ft_fun, w)
    end

    return dft_dw
end