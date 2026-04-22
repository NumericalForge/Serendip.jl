# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


# Solves a system with unknowns in U and F vectors
function solve_system(
    K ::SparseMatrixCSC{Float64, Int},
    U ::Vect,
    F ::Vect,
    nu::Int,
    linear_solver::Symbol=:pardiso
)
    #  ┌  K11   K12 ┐  ┌ U1? ┐    ┌ F1  ┐
    #  │            │  │     │ =  │     │
    #  └  K21   K22 ┘  └ U2  ┘    └ F2? ┘

    linear_solver in (:pardiso, :umfpack) || throw(SerendipException("solve_system: Unknown linear solver $linear_solver"))

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

    # Solve linear system
    F2 = K22*U2
    U1 = zeros(nu)
    if nu>0
        RHS = F1 - K12*U2

        try
            if linear_solver == :pardiso
                ps = MKLPardisoSolver()
                np = get_nprocs(ps)
                nj = Threads.nthreads()
                set_nprocs!(ps, min(np, nj))
                U1 = solve(ps, K11, RHS)
            else
                LUfact = lu(K11)
                U1 = LUfact\RHS
            end

            F2 += K21*U1
        catch err
            err isa InterruptException && rethrow(err)
            if any(isnan.(K11))
                msg = "$msg\nsolve_system: NaN values in coefficients matrix"
            end
            # U1 .= NaN
            return failure("$msg\nsolve_system: $err")
        end
    end

    # maxU = 1e8 # maximum essential value
    maxU = 1/eps() # maximum essential value
    if maximum(abs, U1)>maxU
        return failure("$msg\nsolve_system: Possible syngular matrix ", string(maximum(abs, U1)))
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
    maxits::Int
    rspan::Float64
    alpha::Float64
    beta::Float64
    nmodes::Float64
    eigen_solver::Symbol
    rayleigh::Bool
    linear_solver::Symbol

    @doc """
        SolverSettings(; tol=0.01, rtol=0.01, dT0=0.01, dTmin=1e-7, dTmax=0.1,
                          rspan=0.01, maxits=15, autoinc=false, quiet=false)

    Defines configuration parameters for controlling the analysis process.

    # Arguments
    - `tol::Float64`: Absolute tolerance for convergence checks.
    - `rtol::Float64`: Relative tolerance for convergence checks.
    - `autoinc::Bool`: Enable automatic increment control .
    - `dT0::Float64`: Initial increment of pseudo-time.
    - `dTmin::Float64`: Minimum allowed increment of pseudo-time.
    - `dTmax::Float64`: Maximum allowed increment of pseudo-time.
    - `maxits::Int`: Maximum number of iterations per increment .
    - `rspan::Float64`: Pseudo-time window for reapplying residuals in nonlinear iterations.
    - `alpha::Float64`: Damping coefficient for the mass matrix used in dynamic analyses.
    - `beta::Float64`: Damping coefficient for the stiffness matrix used in dynamic analyses.
    - `nmodes::Int`: Number of modes to compute in modal analysis .
    - `eigen_solver::Symbol`: Modal eigensolver selector (`:auto`, `:arpack`, `:lapack`).
    - `rayleigh::Bool`: Enable Rayleigh damping in dynamic analyses.
    - `linear_solver::Symbol`: Linear solver selector (`:pardiso`, `:umfpack`).
    # Example
    ```julia
    settings = SolverSettings(tol=1e-4, rtol=1e-3, autoinc=true)
    ```
    """
    function SolverSettings(;
        tol=0.01, rtol=0.01, autoinc=false, dT0=0.01, dTmin=1e-7, dTmax=0.1, rspan=0.01,
        maxits=15, alpha=0.0, beta=0.0, nmodes=5, eigen_solver=:auto, rayleigh=false, linear_solver=:pardiso)
        return new(tol, rtol, autoinc, dT0, dTmin, dTmax, maxits, rspan, alpha, beta, nmodes, eigen_solver, rayleigh, linear_solver)
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

    data.log = open("$outdir/analysis.log", "a")

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

        local run_error
        local error_stack
        try
            solstatus = stage_solver(ana, stage, solver_settings; quiet=quiet)
            if succeeded(solstatus)
                stage.status = :done
            else
                stage.status = :failed
            end
        catch err
            run_error = err
            flush(data.log)
            if err isa InterruptException
                stage.status = :interrupted
            else
                stage.status = :error
                error_stack = stacktrace(catch_backtrace())
            end
        end
        close(data.log)

        if !quiet
            wait(status_cycler_task)
            solstatus.message != "" && alert(solstatus.message)
        end

        if stage.status == :interrupted
            throw(SerendipException("The analysis was interrupted"))
        elseif stage.status == :error
            # trim not important frames; try to find the frame that contains REPL/_iterator
            # idx = findfirst(contains("_iterator"), string(frame) for frame in error_stack)
            idx = findfirst(contains("REPL"), string(frame) for frame in error_stack)
            if idx!==nothing
                error_stack = error_stack[1:idx-1]
            end

            alert("Serendip internal error", level=1)
            showerror(stdout, run_error, error_stack)
            # Base.show_backtrace(stdout, error_stack) # shows only the stacktrace
            println()
            stop()
            throw(run_error)
        end

        getlapse(sw)>60 && sound_alert()

    end
    return solstatus

end


"""
    run(ana;
        tol=0.01, rtol=0.01, autoinc=false,
        dT0=0.01, dTmin=1e-7, dTmax=0.1, rspan=0.01, maxits=5,
        alpha=0.0, beta=0.0, nmodes=5, eigen_solver=:auto, rayleigh=false,
        linear_solver=:pardiso,
        quiet=false,
        )

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
- `maxits::Int`: maximum nonlinear iterations per step.
- `alpha::Float64`: Mass matrix coefficient to compute damping in dynamic analyses.
- `beta::Float64`: Stiffness matrix coefficient to compute damping in dynamic analyses.
- `nmodes::Int`: number of modes in modal analysis.
- `eigen_solver::Symbol`: modal eigensolver selector (`:auto`, `:arpack`, `:lapack`).
- `rayleigh::Bool`: enable Rayleigh damping.
- `linear_solver::Symbol`: linear solver selector (`:pardiso`, `:umfpack`).
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
             maxits=10, alpha=0.0, beta=0.25,
             nmodes=8, eigen_solver=:auto, rayleigh=true)
```
"""
function Base.run(ana::Analysis;
    tol     ::Real         = 0.01,
    rtol    ::Real         = 0.01,
    autoinc ::Bool         = false,
    dT0     ::Real         = 0.01,
    dTmin   ::Real         = 1e-7,
    dTmax   ::Real         = 0.2,
    rspan   ::Real         = 0.01,
    maxits  ::Int          = 5,
    alpha   ::Real         = 0.0,
    beta    ::Real         = 0.0,
    nmodes  ::Int          = 5,
    eigen_solver::Symbol     = :auto,
    rayleigh::Bool         = false,
    linear_solver::Symbol  = :pardiso,
    quiet   ::Bool         = false,
    kwargs...,
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
    @check eigen_solver in (:auto, :arpack, :lapack) "run: unknown eigen_solver $eigen_solver"

    if !quiet
        ctx = ana.model.ctx
        printstyled("FE analisys\n", bold=true, color=:cyan)
        println("  type: ", ana.name)
        ctx.stress_state != :auto && println("  stress model: ", ctx.stress_state)
        println("  solver: Newton-Raphson")
        print("  output dir: ")
        printstyled(ana.data.outdir, "\n", color=:cyan)
        print("  output key: ")
        printstyled(ana.data.outkey, "\n", color=:cyan)

        print("  active threads: ")
        nthreads = Threads.nthreads()
        if nthreads==1
            printstyled(Threads.nthreads(), "\n", bold=true, color=:red)
        else
            printstyled(Threads.nthreads(), "\n", color=:green)
        end
    end

    for (key, value) in kwargs
        notify("run: unknown keyword argument $key with value $(repr(value))")
    end

    solver_settings = SolverSettings(
        tol=tol, rtol=rtol, autoinc=autoinc, dT0=dT0, dTmin=dTmin, dTmax=dTmax,
        rspan=rspan, maxits=maxits, alpha=alpha, beta=beta, nmodes=nmodes, 
        eigen_solver=eigen_solver, rayleigh=rayleigh, linear_solver=linear_solver
    )

    status = stage_iterator(ana, solver_settings; quiet=quiet)
    return status
end


function progress_bar(T::Float64)
    dwidth  = displaysize(stdout)[2]
    width   = max(2, min(20, dwidth-21))
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
    infos     = String[]
    alerts    = String[]
    while true
        nlines = 0
        
        _update_msg_list(infos, ana.data.info, 5)
        _update_msg_list(alerts, ana.data.alerts, 5)

        nlines += print_alerts(alerts)
        nlines += print_info(infos)
        nlines += print_summary(ana, sw)
        nlines += print_monitors(ana.data.monitors)

        last_loop && break

        print("\e[$(nlines)A")
        stage.status != :solving && (last_loop=true)
    end

    print("\e[?25h") # enable cursor
end


function _update_msg_list(list::Vector{String}, buff::IOBuffer, maxlen::Int)
    msg = strip(String(take!(buff)))

    if msg!=""
        lines = split(msg, "\n")
        append!(list, lines)
    end

    n = length(list)
    if n>maxlen
        splice!(list, 1:n-maxlen)
        list[1] = "  ⋮"
    end
end


function print_info(infos::Vector{String})
    for m in infos
        printstyled(m, "\e[K\n", color=:cyan)
    end

    return length(infos)
end


function print_alerts(alerts::Vector{String})
    for m in alerts
        printstyled(m, "\e[K\n", color=Base.warn_color())
    end

    return length(alerts)
end


function print_summary(ana::Analysis, sw::StopWatch)
    # display_width  = displaysize(stdout)[2]-2
    data = ana.data 
    nlines = 3

    # line 1:
    T  = data.T
    ΔT = data.ΔT
    printstyled("  inc $(data.inc)  output $(data.out)  udofs $(data.nu)  nf $(data.nfails)\e[K\n", bold=true, color=:light_blue)
    # printstyled("  udofs $(data.nu)", bold=true, color=:light_blue)

    # line 2:
    if data.transient
        t = round(data.t, sigdigits=3)
        printstyled("  t=$t", bold=true, color=:light_blue)
    end
    dT  = round(ΔT,sigdigits=4)
    res_str = data.residue>=0 ? string(round(data.residue,sigdigits=4)) : "…" #"⋯"

    printstyled("  dT $dT  res $res_str\e[K\n", bold=true, color=:light_blue)

    # line 3:
    bar = progress_bar(T)
    progress = @sprintf("%5.3f", T*100)
    printstyled("  $(see(sw)) ", bold=true, color=:light_blue)
    print(bar)
    printstyled(" $(progress)% \e[K\n", bold=true, color=:light_blue)

    return nlines
end


function print_monitors(monitors::Vector{Monitor})
    nlines = 0

    heads  = String[]
    labels = String[]
    values = String[]
    for mon in monitors
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
    for (h, l, v) in zip(heads, labels, values)
        str = "  " * h * " "^(head_len - length(h)) * " " * " "^(label_len - length(l)) * l *" = " * v * "\e[K\n"
        printstyled("  ", str, color=:light_blue)
        nlines += 1
    end

    return nlines
end
    
