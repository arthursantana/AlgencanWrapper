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

#function Dc(i, x)
#    if i == 1
#        return [1, 2], [2*x[1], -1.0]
#    else
#        return [1, 2], [-1.0, -1.0]
#    end
#end
#
#function D²c(i, x)
#    if i == 1
#        return [1], [1], [2]
#    else
#        return [], [], []
#    end
#end

x, lambda = AlgencanWrapper.optimize(
                                     n = 1, m = 2,

                                     f = f, g = ∇f, h = ∇∇f,
                                     equatn = [0, 0],
                                     c = c,# jac = Dc, hc = D²c, equatn = [0, 0],

                                     x = [0.0],
                                     l = [-100.0],
                                     u = [100.0]
                                    )
println(x)
println(f(x))
#println([∇f(x)[1] - lambda[1]*c(1, x)])
#         ∇f(x)[2] - lambda[2]*c(2, x)])

AlgencanWrapper.unload()
