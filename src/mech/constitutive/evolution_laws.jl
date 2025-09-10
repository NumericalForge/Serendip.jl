# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

const _tensile_evolution_laws = ( :constant, :linear, :bilinear, :hordijk, :soft, :custom )
const _compressive_evolution_laws = ( :constant, :default, :popovics )


function setup_compressive_strength(fc::Float64, εc::Float64, fc_law::Union{Symbol,AbstractSpline})
    # fc_law: :constant, :default, :popovics, or custom function (Spline)

    fc < 0.0 || return :none, nothing, failure("fc must be negative. Got fc=$fc)")
    εc < 0.0 || return :none, nothing, failure("εc must be negative. Got εc=$εc)")

    (fc_law in _compressive_evolution_laws || fc_law isa AbstractSpline) ||
        return :none, nothing, failure("fc_law must be :default, :popovics, or a custom function (Spline). Got $(repr(fc_law)).")

    fc_fun = nothing
    if fc_law isa AbstractSpline
        fc_fun = fc_law
        fc_law = :custom

        x0, y0 = fc_fun.X[1], fc_fun.Y[1]
        fc_fun.X[1]!=0.0 && return fc_law, fc_fun, failure("fc_fun must start at (0, fc0). Got ($x0, $y0)")

        # check wether extreme matches fc
        values  = minimum(fc_fun.(range(0.0, fc_fun.X[end], length=40)))
        min_idx = argmin(values)
        fc_     = values[min_idx]
        εc_     = fc_fun.X[min_idx]
        abs((fc_-fc)/fc) < 0.04 || return fc_law, nothing, failure("value of fc from function should match fc value $(fc). Got around $fc_")
        abs(εc_-εc)/abs(εc) < 0.04 || return fc_law, nothing, failure("value of εc from function should match εc value $(εc). Got around $εc_")
    end

    return fc_law, fc_fun, success()
end


function calc_compressive_strength(mat::Constitutive, fc0::Float64, fcr::Float64, εcp::Float64)
    # fc0: initial compressive strength
    # fcr: residual compressive strength
    # εcp: compressive plastic strain

    fc     = mat.fc
    εcp_pk = abs(mat.εc) - abs(fc)/mat.E  # εcp_pk is the compression plastic strain at the peak of the uniaxial compression curve
    χ      = εcp/εcp_pk
    fc_law = mat.fc_law
    # fc0    = 0.35*fc
    # fcr    = 0.1*fc
    
    if fc_law == :constant
        return fc
    elseif fc_law == :popovics
        η = 1/(1 - fc/(mat.E*mat.εc)) # same as  η = 1/(1 - fc/(εcp_pk*Hc))   with   Hc = E*εc/εcp_pk
        return fc * η * χ/(η - 1 + χ^η)
    elseif mat.fc_law==:default
        η = mat.η # from material parameters
        if εcp < 0.0 # sometimes εcp is slightly negative due to numerical errors
            return fc0
        elseif εcp < εcp_pk
            # before peak
            return fc0 + (fc - fc0) * η * χ/(η - 1 + χ^η)
        else
            # after peak
            return fcr + (fc - fcr) * η * χ/(η - 1 + χ^η)
        end
    else
        return mat.fc_fun(εcp)
    end

    return ft
end



function setup_tensile_strength(ft::Float64, GF::Float64, wc::Float64, ft_law::Union{Symbol,AbstractSpline})
    # ft_law: :constant, :hordijk, :linear, :bilinear, :soft, or custom function (Spline)
    # GF: Fracture energy
    # wc: critical crack opening

    ft > 0.0 || return 0.0, :none, nothing, failure("ft must be positive. Got ft=$ft)")

    (ft_law in _tensile_evolution_laws || ft_law isa AbstractSpline) ||
        return wc, ft_law, nothing, failure("ft_law must be :linear, :bilinear, :hordijk, :soft, or a custom function (Spline). Got $(repr(ft_law)).")

    ft_fun = nothing
    if ft_law isa AbstractSpline
        ft_fun = ft_law
        ft_law = :custom

        x0, y0 = ft_fun.X[1], ft_fun.Y[1]
        ft_fun.X[1]!=0.0 && return wc, ft_law, ft_fun, failure("ft_fun must start at (0, tensile_strength=$ft). Got ($x0, $y0)")
        abs(ft_fun.Y[1]-ft)>1e-4 && return wc, ft_law, ft_fun, failure("ft_fun must start at (0, tensile_strength=$ft). Got ($x0, $y0)")

        if ft_fun.Y[end] == 0.0
            wc = ft_fun.X[end]
        else
            wc = Inf
        end
    else
        if isnan(wc)
            isnan(GF) && failure("wc or GF must be defined when using a predefined ft_law model")
            @check GF>0 "GF must be positive. Got GF=$GF"

            if ft_law == :constant
                wc = Inf
            elseif ft_law==:hordijk
                wc = GF/(0.1947*ft)
            elseif ft_law == :linear
                wc = 2*GF/ft
            elseif ft_law == :bilinear
                wc = 5*GF/ft
            elseif ft_law==:soft
                wc = GF/(0.1947*ft)
            end

            wc = round(wc, sigdigits=5)
        end
    end

    wc > 0.0 || return wc, ft_law, nothing, failure("wc must be positive. Got wc=$(repr(wc))")

    return wc, ft_law, ft_fun, success()
end


function calc_tensile_strength(mat::Constitutive, w::Float64)
    # ft: tensile strength at zero crack opening
    # ft_law: :constant, :hordijk, :linear, :bilinear, :soft, or custom function (Spline)
    # ft_fun: custom function (Spline)
    # w: crack opening (plastic displacement jump)

    ft     = mat.ft
    wc     = mat.wc
    ft_law = mat.ft_law

    if ft_law == :constant
        return ft
    elseif ft_law == :hordijk
        if w < wc
            e = exp(1.0)
            z = (1 + 27*(w/wc)^3)*e^(-6.93*w/wc) - 28*(w/wc)*e^(-6.93)           
        else
            z = 0.0
        end
        return z*ft 
    elseif ft_law == :linear
        if w < wc
            a = ft 
            b = ft/wc
        else
            a = 0.0
            b = 0.0
        end
        return a - b*w
    elseif ft_law == :bilinear
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
    elseif ft_law == :soft
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
        @assert mat.ft_fun isa AbstractSpline
        return mat.ft_fun(w)
    end

    return ft
end


function calc_tensile_strength_derivative(mat::Constitutive, w::Float64)
    # ∂σmax/∂up = dft_dw

    ft     = mat.ft
    wc     = mat.wc
    ft_law = mat.ft_law

    if ft_law == :linear
        if w < wc
            b = ft /wc
        else
            b = 0.0
        end
        dft_dw = -b
    elseif ft_law == :bilinear
        σs = 0.25*ft 
        if w < ws
            b  = (ft  - σs)/ws
        elseif w < wc
            b  = σs/(wc-ws)
        else
            b = 0.0
        end
        dft_dw = -b
    elseif ft_law == :hordijk
        if w < wc
            e = exp(1.0)
            dz = ((81*w^2*e^(-6.93*w/wc)/wc^3) - (6.93*(1 + 27*w^3/wc^3)*e^(-6.93*w/wc)/wc) - 0.02738402432/wc)
        else
            dz = 0.0
        end
        dft_dw = dz*ft 
    elseif ft_law == :soft
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
        dft_dw = derive(mat.ft_fun, w)
    end

    return dft_dw
end