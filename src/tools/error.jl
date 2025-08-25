# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export @check


struct StopException <: Exception
    message::String
end

function Base.showerror(io::IO, e::StopException, bt; backtrace=false)
    printstyled(io, e.message, "\n", color=:red)
    # showerror(io, e.message)
    # Base.with_output_color(get(io, :color, false) ? :green : :nothing, io) do io
    #     showerror(io, e.message)
    # end
end

export stop
stop(text="Stop") = throw(StopException(text))


mutable struct SerendipException <: Exception
    message::String
end

Base.showerror(io::IO, e::SerendipException) = printstyled(io, "SerendipException: ", e.message, "\n", color=:red)

# mutable struct PropertyException <: Exception
#     message::String
# end

# Base.showerror(io::IO, e::PropertyException) = print(io, "PropertyException: ", e.message)


macro check(expr, args...)
    exception = ArgumentError
    if length(args)==0
        msg = "Condition `$expr` is not satisfied"
    elseif length(args)==1
        if args[1] isa Exception
            exception = args[1]
        else
            msg = args[1]
        end
    elseif length(args)==2
        exception = args[1]
        msg = args[2]
    else
        error("Invalid number of arguments for @check macro. Expected 0, 1, or 2 arguments, got $(length(args))")
    end

    # expr = replace(string(expr), " "=>"")

    return quote
        if !$(esc(expr)) # Eval boolean expression
            # Get function name
            st = stacktrace(backtrace())
            fname = :_
            for frame in st
                if !startswith(string(frame.func), "_") && frame.func!=Symbol("macro expansion")
                    fname = frame.func
                    break
                end
            end

            # msg = $(esc(msg))
            # msg = repr(begin local value = $(esc(msg)) end)
            msg = $(esc(msg))
            fname != "" && (msg="$fname: $msg")
            throw($(exception)(msg))
        end
    end
end


# macro checkmissing(params, required, names)
#     return quote
#         # Get function name
#         st = stacktrace(backtrace())
#         fname = :_
#         for frame in st
#             if !startswith(string(frame.func), "_") && frame.func!=Symbol("macro expansion")
#                 fname = frame.func
#                 break
#             end
#         end

#         missingkeys = setdiff($(esc(required)), keys($(esc(params))) )
#         if length(missingkeys)>0
#             msg = "Missing arguments: $(join(missingkeys, ", ")). Possible inputs are: $($(esc(names)))"
#             msg = replace(msg, "=" => ":")
#             throw(SerendipException("$fname: $msg"))
#         end
#     end
# end