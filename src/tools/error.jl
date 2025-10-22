# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export @check


struct StopException <: Exception
    message::String
end

function Base.showerror(io::IO, e::StopException, bt; backtrace=false)
    printstyled(io, e.message, "\n", color=:red)
end

export stop
stop(text="Stop") = throw(StopException(text))


mutable struct SerendipException <: Exception
    message::String
end

Base.showerror(io::IO, e::SerendipException) = printstyled(io, "SerendipException: ", e.message, "\n", color=:red)


macro check(expr, args...)
    exception = ArgumentError

    # Default message
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

    ex = Expr(:block)

    # For the case of NaN on parameters
    if expr.head == :call && expr.args[2] isa Symbol && expr.args[1] in (:<, :<=, :>, :>=, :(==), :!=) && expr.args[3] isa Number
        param = expr.args[2]
        msg_nan = "Parameter $param was not provided."
        m = match(r"^(\"?[A-Za-z]\w+)(?=:)", string(msg)) # get the function name if any
        if m !== nothing
            msg_nan = m.match * ": " * msg_nan
        end

        push!(ex.args, :( 
            if isnan($(esc(param))) 
                msg = $(esc(msg_nan))
                throw($(exception)(msg)) 
            end )
        )
    end

    push!(ex.args, :( 
        if !$(esc(expr)) # Eval boolean expression 
            msg = $(esc(msg)) 
            throw($(exception)(msg)) 
        end )
    )

    return ex
end