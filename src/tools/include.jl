# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

include("aliases.jl")

include("error.jl")
include("returnstatus.jl")
include("iteration.jl")

include("show.jl")

include("constants.jl")
include("numerical.jl")
include("splines.jl")
include("linalg.jl")
include("tensors.jl")

include("signal.jl")
include("quaternion.jl")

include("expr.jl") 
include("threads.jl")
include("table.jl")
include("book.jl")
include("utils.jl")
include("stopwatch.jl")
include("xml.jl")
include("encode.jl")


Base.show(io::IO, obj::XmlElement) = _show(io, obj, 3, "")
Base.show(io::IO, obj::XmlDocument)  = _show(io, obj, 3, "")
