using Test

if !isdefined(@__MODULE__, Symbol("@announced_testset"))
    macro announced_testset(args...)
        if length(args) == 2
            name, body = args
            testset_call = Expr(:macrocall, Symbol("@testset"), __source__, name, body)
        elseif length(args) == 3
            opts, name, body = args
            testset_call = Expr(:macrocall, Symbol("@testset"), __source__, opts, name, body)
        else
            error("@announced_testset expects `name begin ... end` or `verbose=true name begin ... end`")
        end

        return esc(quote
            printstyled("\nRunning testset ", string($name), "...\n", color=:light_yellow)
            $testset_call
        end)
    end
end
