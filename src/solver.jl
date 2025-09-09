# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


# Solves a system with unknowns in U and F vectors
function solve_system!(
                       K ::SparseMatrixCSC{Float64, Int},
                       U ::Vect,
                       F ::Vect,
                       nu::Int,
                      )
    #  ┌  K11   K12 ┐  ┌ U1? ┐    ┌ F1  ┐
    #  │            │  │     │ =  │     │
    #  └  K21   K22 ┘  └ U2  ┘    └ F2? ┘

    msg = ""

    # Decomposing the coefficients matrix
    if nu>0
        nu1 = nu+1
        K11 = K[1:nu, 1:nu]
        K12 = K[1:nu, nu1:end]
        K21 = K[nu1:end, 1:nu]
    end
    K22 = K[nu+1:end, nu+1:end]

    F1  = F[1:nu]
    U2  = U[nu+1:end]

    # @showm K11

    # Solve linear system
    F2 = K22*U2
    U1 = zeros(nu)
    # @showm F1
    # @showm F2
    # @showm U2
    if nu>0
        RHS = F1 - K12*U2

        try
            # try
                LUfact = lu(K11)
                U1 = LUfact\RHS
            # catch err
            #     err isa InterruptException && rethrow(err)
            #     if typeof(err)==SingularException
            #         # Regularization attempt
            #         msg = "$msg\nsolve_system!: Syngular matrix - regularization attempt"
            #         S = spdiagm([ 1/maximum(abs, K11[i,:]) for i in 1:nu ])
            #         LUfact = lu(S*K11)
            #         U1  = (LUfact\(S*RHS))
            #     else
            #         return failure("$msg\nsolve_system!: $err")
            #     end
            # end

            F2 += K21*U1
        catch err
            err isa InterruptException && rethrow(err)
            if any(isnan.(K11))
                msg = "$msg\nsolve_system!: NaN values in coefficients matrix"
            end
            # U1 .= NaN
            return failure("$msg\nsolve_system!: $err")
        end
    end

    maxU = 1e8 # maximum essential value
    if maximum(abs, U1)>maxU
        return failure("$msg\nsolve_system!: Possible syngular matrix")
    end

    # Completing vectors
    U[1:nu]     .= U1
    F[nu+1:end] .= F2

    yield()
    return success(msg)
end


struct SolverSettings
    tol::Float64
    rtol::Float64
    autoinc::Bool
    dT0::Float64
    dTmin::Float64
    dTmax::Float64
    tangent_scheme::Symbol
    maxits::Int
    rspan::Float64
    alpha::Float64
    beta::Float64
    nmodes::Float64
    rayleigh::Bool

    @doc """
        SolverSettings(; tol=0.01, rtol=0.01, dT0=0.01, dTmin=1e-7, dTmax=0.1,
                          rspan=0.01, tangent_scheme=:forward_euler, maxits=5, autoinc=false, quiet=false)

    Defines configuration parameters for controlling the analysis process.

    # Arguments
    - `tol::Float64`: Absolute tolerance for convergence checks (default: `0.01`).
    - `rtol::Float64`: Relative tolerance for convergence checks (default: `0.01`).
    - `autoinc::Bool`: Enable automatic increment control (default: `false`).
    - `dT0::Float64`: Initial increment of pseudo-time (default: `0.01`).
    - `dTmin::Float64`: Minimum allowed increment of pseudo-time (default: `1e-7`).
    - `dTmax::Float64`: Maximum allowed increment of pseudo-time (default: `0.1`).
    - `tangent_scheme::Symbol`: Tangent update approach (`:forward_euler`, `:heun`, `ralston`).
    - `maxits::Int`: Maximum number of iterations per increment (default: `5`).
    - `rspan::Float64`: Progression span for reapplying the residual in nonlinear iterations (default: `0.01`).
    - `alpha::Float64`: Damping coefficient for the mass matrix used in dynamic analyses (default: `0.0`).
    - `beta::Float64`: Damping coefficient for the stiffness matrix used in dynamic analyses (default: `0.0`).
    - `nmodes::Int`: Number of modes to compute in modal analysis (default: `5`).

    # Example
    ```julia
    settings = SolverSettings(tol=1e-4, rtol=1e-3, autoinc=true)
    ```
    """
    function SolverSettings(;
        tol=0.01, rtol=0.01, autoinc=false, dT0=0.01, dTmin=1e-7, dTmax=0.1, rspan=0.01,
        tangent_scheme=:forward_euler, maxits=5, alpha=0.0, beta=0.0, nmodes=5, rayleigh=false)
        return new(tol, rtol, autoinc, dT0, dTmin, dTmax, tangent_scheme, maxits, rspan, alpha, beta, nmodes, rayleigh)
    end
end


function stage_iterator(ana::Analysis, solver_settings::SolverSettings; quiet::Bool=false)
    autoinc = solver_settings.autoinc
    data    = ana.data

    cstage = findfirst(st->st.status!=:done, ana.data.stages)
    cstage === nothing && throw(SerendipException("stage_iterator: No stages have been set for $(ana.name)"))

    solstatus = success()

    outdir = ana.data.outdir

    if !isdir(outdir)
        info("solve!: creating output directory ./$outdir")
        mkpath(outdir)
    end


    if cstage==1
        data.log = open("$outdir/solve.log", "w")
    else
        data.log = open("$outdir/solve.log", "a")
    end


    for stage in ana.data.stages[cstage:end]
        stage.status = :solving

        nincs  = stage.nincs
        nouts  = stage.nouts

        data.stage = stage.id
        data.inc   = 0
        data.T = 0.0

        if !quiet
            printstyled("Stage $(stage.id)\n", bold=true, color=:cyan)
        end

        save_outs = stage.nouts > 0
        if save_outs
            if !autoinc
                if nouts > nincs
                    nincs = nouts
                    quiet || info("nincs changed to $(nincs) to match nouts")
                end
                if nincs%nouts != 0
                    stage.nincs = nincs - (nincs%nouts) + nouts
                    quiet || info("nincs changed to $nincs to be multiple of nouts")
                end
            end
            stage.nincs = nincs
            stage.nouts = nouts
        end

        sw = StopWatch() # timing
        if !quiet
            status_cycler_task = Threads.@spawn :interactive status_cycler(ana, sw)
        end

        local runerror
        local error_st
        try
            solstatus = stage_solver(ana, stage, solver_settings; quiet=quiet)
            if succeeded(solstatus)
                stage.status = :done
            else
                stage.status = :failed
            end
        catch err
            runerror = err
            flush(data.log)
            if err isa InterruptException
                stage.status = :interrupted
            else
                stage.status = :error
                error_st = stacktrace(catch_backtrace())
            end
        end
        close(data.log)

        if !quiet
            wait(status_cycler_task)
            solstatus.message != "" && println(solstatus.message)
        end

        if stage.status == :interrupted
            throw(SerendipException("The analysis was interrupted"))
        elseif stage.status == :error
            # trim not important frames; try to find the frame that contains REPL/_iterator
            # idx = findfirst(contains("_iterator"), string(frame) for frame in error_st)
            idx = findfirst(contains("REPL"), string(frame) for frame in error_st)
            if idx!==nothing
                error_st = error_st[1:idx-1]
            end

            alert("Serendip internal error", level=1)
            showerror(stdout, runerror, error_st)
            # Base.show_backtrace(stdout, error_st) # shows only the stacktrace
            println()
            stop()
            throw(runerror)
        end

        getlapse(sw)>60 && sound_alert()

    end
    return solstatus

end


"""
    run(ana;
        tol=0.01, rtol=0.01, autoinc=false,
        dT0=0.01, dTmin=1e-7, dTmax=0.1, rspan=0.01,
        tangent_scheme=:forward_euler, maxits=5,
        alpha=0.0, beta=0.0, nmodes=5, rayleigh=false,
        quiet=false)

Execute a finite-element analysis and return the solver status.

# Arguments
- `ana::Analysis`: analysis object to solve.

# Keywords
- `tol::Real`: absolute convergence tolerance for the residual.
- `rtol::Real`: relative convergence tolerance for the residual.
- `autoinc::Bool`: enable automatic step size control.
- `dT0::Real`: initial time/load increment.
- `dTmin::Real`: minimum allowed increment.
- `dTmax::Real`: maximum allowed increment.
- `rspan::Real`: span parameter for the auto-increment controller.
- `tangent_scheme::Symbol`: global tangent computation approach (`:forward_euler`, `:heun`, `:ralston`).
- `maxits::Int`: maximum nonlinear iterations per step.
- `alpha::Float64`: Mass matrix coefficient to compute damping in dynamic analyses.
- `beta::Float64`: Stiffness matrix coefficient to compute damping in dynamic analyses.
- `nmodes::Int`: number of modes in modal analysis.
- `rayleigh::Bool`: enable Rayleigh damping.
- `quiet::Bool`: suppress console output.

# Behavior
- Prints a short banner with analysis type, stress model, scheme, and active threads unless `quiet=true`.
- Builds `SolverSettings` from the provided keywords and advances stages via `stage_iterator`.

# Returns
- Solver status object returned by `stage_iterator` (implementation-specific).

# Example
```julia
status = run(analysis;
             tol=1e-3, rtol=1e-3, autoinc=true,
             dT0=0.02, dTmin=1e-6, dTmax=0.1,
             tangent_scheme=:rk2, maxits=10, alpha=0.0, beta=0.25,
             nmodes=8, rayleigh=true)
```
"""
function Base.run(ana::Analysis;
    tol     ::Real         = 0.01,
    rtol    ::Real         = 0.01,
    autoinc ::Bool         = false,
    dT0     ::Real         = 0.01,
    dTmin   ::Real         = 1e-7,
    dTmax   ::Real         = 0.1,
    rspan   ::Real         = 0.01,
    scheme  ::Symbol       = :none,
    tangent_scheme::Symbol = :forward_euler,
    maxits  ::Int          = 5,
    alpha   ::Real         = 0.0,
    beta    ::Real         = 0.0,
    nmodes  ::Int          = 5,
    rayleigh::Bool         = false,
    quiet   ::Bool         =false
)
    @check tol>0 "run: solver paramter `tol` must be positive"
    @check rtol>0 "run: solver paramter `rtol` must be positive"
    @check dT0>0 "run: solver paramter `dT0` must be positive"
    @check dTmin>0 "run: solver paramter `dTmin` must be positive"
    @check dTmax>0 "run: solver paramter `dTmax` must be positive"
    @check dTmin<dTmax "run: solver paramter `dTmin` must be less than `dTmax`"
    @check rspan>0 "run: solver paramter `rspan` must be positive"
    @check maxits>0 "run: solver paramter `maxits` must be positive"
    @check nmodes>0 "run: solver paramter `nmodes` must be positive"
    @check alpha>=0 "run: solver paramter `alpha` must be positive"
    @check beta>=0 "run: solver paramter `beta` must be positive"

    if scheme != :none # for backward compatibility
        alert("run: 'scheme' keyword is deprecated; use 'tangent_scheme' with options :forward_euler, :heun, :ralston")
        tangent_scheme == scheme==:ME ? :heun : scheme == :FE ? :forward_euler : tangent_scheme
    end
    @check tangent_scheme in (:forward_euler, :heun, :ralston) "run: unknown tangent_scheme $tangent_scheme"

    if !quiet
        ctx = ana.model.ctx
        printstyled("FE analisys\n", bold=true, color=:cyan)
        println("  type: ", ana.name)
        ctx.stress_state != :auto && println("  stress model: ", ctx.stress_state)
        println("  tangent rule: ", tangent_scheme)

        print("  active threads: ")
        nthreads = Threads.nthreads()
        if nthreads==1
            printstyled(Threads.nthreads(), "\n", color=:red)
        else
            printstyled(Threads.nthreads(), "\n", color=:green)
        end
    end


    solver_settings = SolverSettings(
        tol=tol, rtol=rtol, autoinc=autoinc, dT0=dT0, dTmin=dTmin, dTmax=dTmax,
        rspan=rspan, tangent_scheme=tangent_scheme, maxits=maxits, alpha=alpha, beta=beta, nmodes=nmodes, rayleigh=rayleigh,
    )

    status = stage_iterator(ana, solver_settings; quiet=quiet)
    return status
end


function progress_bar(T::Float64)
    dwidth  = displaysize(stdout)[2]-2
    width   = max(2, min(20, dwidth-35))
    ch_done = T*width
    frac    = ch_done - floor(ch_done)

    barl = repeat(['━'], floor(Int, ch_done))
    barr = Char[]

    if frac<=0.25
        push!(barr, '╶')
        push!(barr, '─')
    elseif 0.25<frac<0.75
        push!(barl, '╸')
        push!(barr, '─')
    else
        push!(barl, '━')
        push!(barr, '╶')
    end

    if length(barl)>=width
        barl = barl[1:width]
        barr = Char[]
    end

    if length(barl)+length(barr)==width+1
        barr = Char[barr[1]]
    end

    append!(barr, repeat(['─'], width -length(barl) -length(barr) ))
    barls = reduce(*, barl)
    barrs = reduce(*, barr)

    iscolor = get(stdout, :color, false)
    if iscolor
        color        = :blue
        enable_color = get(Base.text_colors, color, Base.text_colors[:default])
        enable_bold  = get(Base.text_colors, :bold, Base.text_colors[:default])
        normal_color = get(Base.disable_text_style, :normal, Base.text_colors[:default])
        disable_bold = get(Base.disable_text_style, :bold, Base.text_colors[:default])
        barls        = string(enable_color, enable_bold, barls, disable_bold, normal_color)
        enable_color = get(Base.text_colors, :light_black, Base.text_colors[:default])
        normal_color = get(Base.disable_text_style, :bold, Base.text_colors[:default])
        barrs        = string(enable_color, barrs, normal_color)
    end

    return barls*barrs

end


function status_cycler(ana::Analysis, sw::StopWatch)
    print("\e[?25l") # disable cursor

    stage     = ana.data.stages[ana.data.stage]
    last_loop = false
    alerts    = String[]
    while true
        nlines = 0

        nlines += print_info(ana)
        nlines += print_alerts(ana, alerts)
        nlines += print_summary(ana, sw)

        last_loop && break

        print("\e[$(nlines)A")
        stage.status != :solving && (last_loop=true)
        sleep(0.05)
        # yield()
    end

    print("\e[?25h") # enable cursor
end


function print_info(ana::Analysis)
    str = strip(String(take!(ana.data.info)))
    str!="" && println("  ", str, "\e[K")

    return 0
end


function print_alerts(ana::Analysis, alerts::Array{String,1})
    str = strip(String(take!(ana.data.alerts)))

    if str!=""
        list = split(str, "\n")
        list = String[ string("  ", Time(now()), "  ", m) for m in list ]
        append!(alerts, list)
    end

    n = length(alerts)
    if n>5
        splice!(alerts, 1:n-5)
        alerts[1] = "  ⋮"
    end

    for m in alerts
        printstyled(m, "\e[K\n", color=Base.warn_color())
    end

    return length(alerts)
end


function print_summary(ana::Analysis, sw::StopWatch)
    data = ana.data
    nlines = 2

    # line 1:
    T  = data.T
    ΔT = data.ΔT
    printstyled("  inc $(data.inc) output $(data.out)", bold=true, color=:light_blue)
    if data.transient
        t = round(data.t, sigdigits=3)
        printstyled(" t=$t", bold=true, color=:light_blue)
    end
    dT  = round(ΔT,sigdigits=4)
    res = round(data.residue,sigdigits=4)

    printstyled(" dT=$dT res=$res\e[K\n", bold=true, color=:light_blue)

    # line 2:
    bar = progress_bar(T)
    progress = @sprintf("%5.3f", T*100)
    printstyled("  $(see(sw)) ", bold=true, color=:light_blue)
    print(bar)
    printstyled(" $(progress)% \e[K\n", bold=true, color=:light_blue)

    # print monitors
    heads  = String[]
    labels = String[]
    values = String[]
    for mon in ana.data.monitors
        h, ls, vs = output(mon)
        length(ls) > 0 || continue
        hs = repeat([""], length(ls))
        hs[1] = h * " :"
        append!(heads, hs)
        append!(labels, ls)
        append!(values, vs)
    end

    head_len = maximum(length.(heads), init=0)
    label_len = maximum(length.(labels), init=0)
    print("\e[K") # clear line
    for (h, l, v) in zip(heads, labels, values)
        str = "  " * h * " "^(head_len - length(h)) * " " * " "^(label_len - length(l)) * l *" = " * v * "\e[K\n"
        printstyled("  ", str, color=:light_blue)
        nlines += 1
    end

    return nlines
end
