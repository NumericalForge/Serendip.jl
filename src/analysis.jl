# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

"""
    AnalysisData(; outkey="", outdir="")

Container for global data produced and updated during a numerical analysis.
Holds information about stages, outputs, logging, monitors, time control, and solver state.

# Fields
- `stages::Vector{Stage}`: List of analysis stages.
- `loggers::Vector{Logger}`: Active loggers for outputs.
- `monitors::Vector{Monitor}`: Active monitors for convergence, criteria, or variables.

- `outdir::String`: Output directory (default = `"./output"`).
- `outkey::String`: Key prefix for output files (default = `"out"`).

- `log::IOStream`: Stream for solver log file (`analysis.log` in `outdir`).
- `alerts::IOBuffer`: Buffer for warnings and alerts.
- `info::IOBuffer`: Buffer for informational messages.

- `T::Float64`: Current pseudo–time for stage progression.
- `t::Float64`: Current physical time (for transient analysis).
- `transient::Bool`: Flag for time‐dependent analysis.

- `ΔT::Float64`: Pseudo–time increment of the current step.
- `Tupdate::Float64`: Pseudo–time of last output update.
- `flushtime::Float64`: Time of last flush to output.

- `stage::Int`: Current stage counter.
- `inc::Int`: Current increment counter.
- `out::Int`: Current output file counter.

- `residue::Float64`: Convergence residue for the current iteration.

# Constructor
```julia
AnalysisData(; outkey="", outdir="")
```
Creates a new `AnalysisData` instance with specified output key and directory.
If `outkey` is empty, defaults to `"out"`. If `outdir` is empty, defaults to `"./output"`.
The output directory is created if it does not exist.
"""
mutable struct AnalysisData
    stages::Vector{Stage}
    loggers::Vector{Logger}
    monitors::Vector{Monitor}

    outdir   ::String    # Output directory
    outkey   ::String    # Key name for output files

    log      ::IOStream # solver log file
    alerts   ::IOBuffer # alerts
    info     ::IOBuffer

    T      ::Float64  # Pseudo time for current stage
    t      ::Float64   # Time in time dependent analysis
    transient::Bool      # Time dependent analysis

    ΔT     ::Float64  # Pseudo time current increment
    Tupdate::Float64  # Pseudo time for last files update
    flushtime::Float64 # Time of last flush

    # Current counters
    stage  ::Int      # Current stage
    inc    ::Int      # Current increment
    out    ::Int      # Current output file number
    nu     ::Int      # Current number of unknown dofs

    residue::Float64  # Residue

    function AnalysisData(; outkey="", outdir="")
        this = new()
        outkey = isempty(outkey) ? "out" : outkey
        outdir = isempty(outdir) ? "./output" : rstrip(outdir, ['/', '\\'])
        isdir(outdir) || mkpath(outdir)

        this.stages    = Stage[]
        this.loggers   = Logger[]
        this.monitors  = Monitor[]
        this.outkey    = outkey
        this.outdir    = outdir
        this.log       = open(joinpath(outdir, "analysis.log"), "w")
        this.alerts    = IOBuffer()
        this.info      = IOBuffer()
        this.T         = 0.0
        this.t         = 0.0
        this.transient = false
        this.ΔT        = 0.0
        this.Tupdate   = 0.0
        this.flushtime = 0.0
        this.stage     = 0
        this.inc       = 0
        this.out       = 0
        this.nu        = 0
        this.residue   = 0.0

        return this
    end
end


function add_stage(
    ana::Analysis;
    name::String="",
    nincs::Int=1,
    nouts::Int=0,
    tspan::Float64=0.0,
    activate::Any=:all,  # filtert to activate
    deactivate::Any=:all # filtert to deactivate
    )

    _activate   = select(ana.model, :element, activate)
    _deactivate = select(ana.model, :element, deactivate)

    stage = Stage(name; nincs=nincs, nouts=nouts, tspan=tspan, activate=_activate, deactivate=_deactivate)
    stage.id = length(ana.data.stages) + 1
    stage.analysis = ana
    push!(ana.data.stages, stage)
    return stage
end


function compute_bc_values(ana::Analysis, bc::BoundaryCondition, t::Float64, U::Vector{Float64}, F::Vector{Float64})

    if bc.kind == :node
        # essential_keys = Set( dof.name for node in bc.nodes for dof in node.dofs )
        # ess_keys = get_essential_keys(ana)
        # nat_keys = get_natural_keys(ana)
        for node in bc.target
            ess_keys = get_essential_keys(node)
            nat_keys = get_natural_keys(node)
            x, y, z = node.coord
            for (key,cond) in bc.conds
                if key in ess_keys
                    dof = get_dof(node, key)
                    U[dof.eq_id] = evaluate(cond, x=x, y=y, z=z, t=t)
                elseif key in nat_keys
                    dof = get_dof(node, key)
                    F[dof.eq_id] += evaluate(cond, x=x, y=y, z=z, t=t) # cummulative value
                else
                    @warn("compute_bc_values: Unknown boundary condition key `$(repr(key))`. Available keys: $(ess_keys ∪ nat_keys).")
                end
            end
        end
    elseif bc.kind in (:face, :edge)
        facets = bc.target
        essential_keys = Set( dof.name for facet in facets for node in facet.nodes for dof in node.dofs )

        for (key,cond) in bc.conds
            if key in essential_keys
                for facet in facets
                    for node in facet.nodes
                        dof = get_dof(node, key)
                        dof===nothing && continue
                        x, y, z = node.coord
                        U[dof.eq_id] = evaluate(cond, x=x, y=y, z=z, t=t)
                    end
                end
            else
                for facet in facets
                    Fd, map = distributed_bc(facet.owner, facet, t, key, cond)
                    F[map] += Fd
                end
            end
        end
    elseif bc.kind == :body
        for elem in bc.target
            for (key,val) in bc.conds
                Fd, map = body_load(elem, key, val)
                F[map] += Fd
            end
        end
    else
        error("Unsupported boundary condition kind: $(bc.kind). Available kinds: :node, :face, :edge, :body")
    end

end


# Return a vector with all model dofs and the number of unknown dofs according to bcs
function configure_dofs(model::AbstractDomain, bcs::Vector{BoundaryCondition})

    # get active nodes
    # active_elems = select(model.elems, :active)
    # ids = [ node.id for elem in active_elems for node in elem.nodes ]
    # ids = sort(unique(ids)) # sort is required to preserve node numbering optimization

    ids = [ node.id for node in model.nodes ]

    active_nodes = model.nodes[ids]

    # All dofs
    # dofs = Dof[dof for node in active_nodes for dof in node.dofs]
    dofs = Dof[dof for node in active_nodes if !node.aux for dof in node.dofs]

    # Reset all dofs as natural conditions
    for dof in dofs
        dof.prescribed = false
    end

    for bc in bcs
        configure_bc_dofs(bc)
    end

    # Split dofs
    presc = [ dof.prescribed for dof in dofs ]
    pdofs = dofs[presc]
    udofs = dofs[.!presc]
    dofs  = [ udofs; pdofs ]
    nu    = length(udofs)

    # set eq_id in dofs
    for (i,dof) in enumerate(dofs)
        dof.eq_id = i
    end

    return dofs, nu
end


function update_records!(ana::Analysis; checkpoint=true, force=false)
    data = ana.data
    outdir = data.outdir

    flushinterval = 5.0
    flush = time()-data.flushtime>flushinterval || checkpoint || force || data.T >= 1.0-1e-8
    flush && (data.flushtime = time())

    if checkpoint
        rm.(glob("*conflicted*.log"), force=true)
        rm.(glob("*conflicted*.*", "$outdir/"), force=true)

        update_output_data!(ana.model) # need to be before group loggers
        save(ana.model, "$outdir/$(data.outkey)-$(data.out).vtu", quiet=true)

        # update multiloggers
        for logger in ana.data.loggers
            if logger.kind in (:ipgroup, :nodegroup)
                update_logger!(logger, ana)
                logger.filename!="" && save(logger.table, logger.filename, quiet=true)
            end
        end

    end

    flush && Base.flush(data.log)

    # update single loggers
    data.Tupdate > data.T && (data.Tupdate=0.0) # for subsequent stages
    update_single_loggers = data.T-data.Tupdate >= 0.00025 || data.T==0
    update_single_loggers && (data.Tupdate = data.T)

    for logger in ana.data.loggers
        if logger.kind in (:node, :ip, :nodalreduce)
            update_single_loggers && update_logger!(logger, ana)
            flush && logger.filename!="" && save(logger.table, logger.filename, quiet=true)
        end
    end

    # update monitors
    for monitor in ana.data.monitors
        rstatus = update_monitor!(ana, monitor)
        failed(rstatus) && return rstatus
        flush && monitor.filename!="" && save(monitor.table, monitor.filename, quiet=true)
    end

    return success()
end