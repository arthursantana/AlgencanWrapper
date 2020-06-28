include("AlgencanWrapper.jl")

function f(x)
    return -(x[1]-1)^2
end

function ∇f(x)
    return [-2x[1] + 2]
end

function ∇∇f(x)
    return [1], [1], [-2]
end

function c(i, x)
    if i == 1
        return -x[1] - 10
    else
        return x[1] + 10
    end
end

x, lambda = AlgencanWrapper.optimize(
                                     n = 1, m = 2,

                                     f = f, g = ∇f, h = ∇∇f,
                                     equatn = [0, 0],
                                     c = c,

                                     x = [0.0],
                                     l = [-100.0],
                                     u = [100.0],
                                    )
println(x)
println(f(x))

AlgencanWrapper.unload()
