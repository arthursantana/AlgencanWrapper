include("AlgencanWrapper.jl")

function f(x)
    return x[2]
end

function ∇f(x)
    return [0, 1]
end

function ∇∇f(x)
    return [], [], []
end

function c(i, x)
    if i == 1
        return x[1]^2 + 1.0 - x[2]
    else
        return 2.0 - x[1] - x[2]
    end
end

function Dc(i, x)
    if i == 1
        return [1, 2], [2*x[1], -1.0]
    else
        return [1, 2], [-1.0, -1.0]
    end
end

function D²c(i, x)
    if i == 1
        return [1], [1], [2]
    else
        return [], [], []
    end
end

x, fx = AlgencanWrapper.optimize(2, f, ∇f, ∇∇f, 2, c, [0 0], Dc, D²c)
println(x)
println(f(x))

AlgencanWrapper.unload()
