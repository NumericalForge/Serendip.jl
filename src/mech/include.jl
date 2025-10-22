# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


include("elements/mech.jl")
include("elements/utils.jl")
include("elements/distributed.jl")

# Elements
include("elements/mech-frame.jl")
include("elements/mech-bulk.jl")
# include("elements/mech-fluid.jl")
include("elements/mech-bar.jl")
include("elements/mech-embtruss.jl")
include("elements/mech-beam.jl")
include("elements/mech-contact.jl")
include("elements/mech-cohesive.jl")
include("elements/mech-bondslip.jl")
include("elements/mech-bondtip.jl")
include("elements/mech-shell.jl")

include("constitutive/evolution_laws.jl")

# Models for bulk elements
include("constitutive/linear-elastic.jl")
# include("constitutive/linear-elastic-fluid.jl")
include("constitutive/drucker-prager.jl")
include("constitutive/von-mises.jl")
include("constitutive/willam-warnke.jl")
include("constitutive/UCP.jl")

# Spring, dumper, lumped mass
include("elements/mech-lumpedmass.jl")
include("constitutive/lumpedmass.jl")
include("elements/mech-spring.jl")
include("constitutive/linear-spring.jl")

# Models for interface, interface and coohesive elements
include("constitutive/elastic-interface.jl") # includes: LinearIntertace, LinearContact
include("constitutive/mohr-coulomb-contact.jl")
include("constitutive/linear-cohesive.jl")
include("constitutive/mohr-coulomb-cohesive.jl")
include("constitutive/power-yield-cohesive.jl")
include("constitutive/asinh-yield-cohesive.jl")

# Models for 1D interface elements
include("constitutive/linear-bondslip.jl")
include("constitutive/ceb-bondslip.jl")
include("constitutive/power-exp-bondslip.jl")
include("constitutive/cyclic-bondslip.jl")
include("constitutive/elastic-tip.jl") # includes: LinearTip, LinearTipContact

# Stress-update integrator
include("constitutive/integrator.jl")

# Solvers
include("mech-analysis.jl")
include("dyn-analysis.jl")
include("mech-modal-analysis.jl")
