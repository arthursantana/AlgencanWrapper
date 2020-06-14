module AlgencanWrapper

using Libdl

algencan = Libdl.dlopen("./libalgencan")
algencan_optimize = Libdl.dlsym(algencan, "c_algencan")

function unsafe_store_array!(dest, convert_to, src)
    for (i, v) in enumerate(src)
        unsafe_store!(dest, convert_to(v), i)
    end
end

function make_evalf(fun)
    return function (n::Cint, x::Ptr{Cdouble}, f::Ptr{Cdouble}, flag::Ptr{Cint})
        if fun == nothing
            unsafe_store!(flag, Cint(-1))
            return nothing
        end

        j_x = unsafe_wrap(Array, x, Int(n))

        unsafe_store!(f, Cdouble(fun(j_x)))

        unsafe_store!(flag, Cint(0))
        return nothing
    end
end

function make_evalg(fun)
    return function (n::Cint, x::Ptr{Cdouble}, g::Ptr{Cdouble}, flag::Ptr{Cint})
        if fun == nothing
            unsafe_store!(flag, Cint(-1))
            return nothing
        end

        j_x = unsafe_wrap(Array, x, Int(n))

        unsafe_store_array!(g, Cdouble, fun(j_x))

        unsafe_store!(flag, Cint(0))
        return nothing
    end
end

function make_evalh(fun)
    return function (n::Cint, x::Ptr{Cdouble},
                     hrow::Ptr{Cint}, hcol::Ptr{Cint},
                     hval::Ptr{Cdouble}, hnnz::Ptr{Cint},
                     lim::Cint, lmem::Ptr{Cuchar}, flag::Ptr{Cint})
        if fun == nothing
            unsafe_store!(flag, Cint(-1))
            return nothing
        end

        j_x = unsafe_wrap(Array, x, Int(n))

        j_hrow, j_hcol, j_hval = fun(j_x)
        j_hnnz = length(j_hval)

        if j_hnnz > Int(lim)
            unsafe_store!(lmem, Cuchar(1))
            unsafe_store!(flag, Cint(-1))
            return nothing
        end

        unsafe_store_array!(hrow, Cint, j_hrow)
        unsafe_store_array!(hcol, Cint, j_hcol)
        unsafe_store_array!(hval, Cdouble, j_hval)
        unsafe_store!(hnnz, Cint(j_hnnz))

        unsafe_store!(lmem, Cuchar(0))

        unsafe_store!(flag, Cint(0))
        return nothing
    end
end

function make_evalc(fun)
    return function (n::Cint, x::Ptr{Cdouble}, ind::Cint, c::Ptr{Cdouble}, flag::Ptr{Cint})
        if fun == nothing
            unsafe_store!(flag, Cint(-1))
            return nothing
        end

        j_n = Int(n)
        j_ind = Int(ind) + 1 # convert 0-indexed to 1-indexed

        if j_ind > j_n
            unsafe_store!(flag, Cint(-1))
            return nothing
        end

        j_x = unsafe_wrap(Array, x, j_n)

        j_c = fun(j_ind, j_x)

        unsafe_store!(c, Cdouble(j_c))
        unsafe_store!(flag, Cint(0))
        return nothing
    end
end


function make_evaljac(fun)
    return function (n::Cint, x::Ptr{Cdouble}, ind::Cint,
                     jcvar::Ptr{Cint}, jcval::Ptr{Cdouble}, jcnnz::Ptr{Cint},
                     lim::Cint, lmem::Ptr{Cuchar}, flag::Ptr{Cint})
        if fun == nothing
            unsafe_store!(flag, Cint(-1))
            return nothing
        end

        j_n = Int(n)
        j_ind = Int(ind) + 1 # convert 0-indexed to 1-indexed

        if j_ind > j_n
            unsafe_store!(flag, Cint(-1))
            return nothing
        end

        j_x = unsafe_wrap(Array, x, Int(n))

        j_jcvar, j_jcval = fun(j_ind, j_x)
        j_jcnnz = length(j_jcval)

        if j_jcnnz > Int(lim)
            unsafe_store!(lmem, Cuchar(1))
            unsafe_store!(flag, Cint(-1))
            return nothing
        end

        unsafe_store_array!(jcvar, Cint, map(x -> x - 1, j_jcvar)) # convert back
        unsafe_store_array!(jcval, Cdouble, j_jcval)
        unsafe_store!(jcnnz, Cint(j_jcnnz))

        unsafe_store!(lmem, Cuchar(0))
        unsafe_store!(flag, Cint(0))
        return nothing
    end
end

function make_evalhc(fun)
    return function (n::Cint, x::Ptr{Cdouble}, ind::Cint,
                     hcrow::Ptr{Cint}, hccol::Ptr{Cint},
                     hcval::Ptr{Cdouble}, hcnnz::Ptr{Cint},
                     lim::Cint, lmem::Ptr{Cuchar}, flag::Ptr{Cint})
        if fun == nothing
            unsafe_store!(flag, Cint(-1))
            return nothing
        end

        j_n = Int(n)
        j_ind = Int(ind) + 1 # convert 0-indexed to 1-indexed

        if j_ind > j_n
            unsafe_store!(flag, Cint(-1))
            return nothing
        end

        j_x = unsafe_wrap(Array, x, Int(n))

        j_hcrow, j_hccol, j_hcval = fun(j_ind, j_x)
        j_hcnnz = length(j_hcval)

        if j_hcnnz > Int(lim)
            unsafe_store!(lmem, Cuchar(1))
            unsafe_store!(flag, Cint(-1))
            return nothing
        end

        unsafe_store_array!(hcrow, Cint, map(x -> x - 1, j_hcrow)) # convert back
        unsafe_store_array!(hccol, Cint, map(x -> x - 1, j_hccol)) # convert back
        unsafe_store_array!(hcval, Cdouble, j_hcval)
        unsafe_store!(hcnnz, Cint(j_hcnnz))

        unsafe_store!(lmem, Cuchar(0))
        unsafe_store!(flag, Cint(0))
        return nothing
    end
end

function myevalfc(n::Cint, x::Ptr{Cdouble}, f::Ptr{Cdouble}, m::Cint, c::Ptr{Cdouble}, flag::Ptr{Cint})
    unsafe_store!(flag, Cint(-1))
    return nothing
end

function myevalgjac(n::Cint, x::Ptr{Cdouble}, g::Ptr{Cdouble}, m::Cint, jcfun::Ptr{Cint}, jcvar::Ptr{Cint},
                    jcval::Ptr{Cdouble}, jcnnz::Ptr{Cint}, lim::Cint, lmem::Ptr{Cuchar}, flag::Ptr{Cint})
    unsafe_store!(flag, Cint(-1))
    return nothing
end

function myevalgjacp(n::Cint, x::Ptr{Cdouble}, g::Ptr{Cdouble}, m::Cint, p::Ptr{Cdouble}, q::Ptr{Cdouble},
                     work::Cchar, gotj::Ptr{Cuchar}, flag::Ptr{Cint})
    unsafe_store!(flag, Cint(-1))
    return nothing
end

function myevalhl(n::Cint, x::Ptr{Cdouble}, m::Cint, lambda::Ptr{Cdouble}, scalef::Cdouble,
                  scalec::Ptr{Cdouble}, hlrow::Ptr{Cint}, hlcol::Ptr{Cint}, hlval::Ptr{Cdouble},
                  hlnnz::Ptr{Cint}, lim::Cint, lmem::Ptr{Cuchar}, flag::Ptr{Cint})
    unsafe_store!(flag, Cint(-1))
    return nothing
end

function myevalhlp(n::Cint, x::Ptr{Cdouble}, m::Cint, lambda::Ptr{Cdouble}, scalef::Cdouble,
                   scalec::Ptr{Cdouble}, p::Ptr{Cdouble}, hp::Ptr{Cdouble},
                   goth::Ptr{Cuchar}, flag::Ptr{Cint})
    unsafe_store!(flag, Cint(-1))
    return nothing
end

function optimize(n::Int,
                  f::Function,
                  g::Union{Function, Nothing}=nothing,
                  h::Union{Function, Nothing}=nothing,
                  m::Int=0,
                  c::Union{Function, Nothing}=nothing,
                  equatn::Union{Array{Int}, Nothing}=nothing,
                  jac::Union{Function, Nothing}=nothing,
                  hc::Union{Function, Nothing}=nothing)
    myevalf = make_evalf(f)
    myevalg = make_evalg(g)
    myevalh = make_evalh(h)
    myevalc = make_evalc(c)
    myevaljac = make_evaljac(jac)
    myevalhc = make_evalhc(hc)

    myevalf_c = @cfunction($myevalf, Cvoid, (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}))
    myevalg_c = @cfunction($myevalg, Cvoid, (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}))
    myevalh_c = @cfunction($myevalh, Cvoid, (Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
                                             Ptr{Cint}, Cint, Ptr{Cuchar}, Ptr{Cint}))
    myevalc_c = @cfunction($myevalc, Cvoid, (Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cint}))
    myevaljac_c = @cfunction($myevaljac, Cvoid, (Cint, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cdouble},
                                                 Ptr{Cint}, Cint, Ptr{Cuchar}, Ptr{Cint}))
    myevalhc_c = @cfunction($myevalhc, Cvoid, (Cint, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Cint, Ptr{Cuchar}, Ptr{Cint}))
    myevalfc_c = @cfunction($myevalfc, Cvoid, (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cint}))
    myevalgjac_c = @cfunction($myevalgjac, Cvoid, (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint},
                                                   Ptr{Cint},
                                                   Ptr{Cdouble}, Ptr{Cint}, Cint, Ptr{Cuchar}, Ptr{Cint}))
    myevalgjacp_c = @cfunction($myevalgjacp, Cvoid, (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble},
                                                     Ptr{Cdouble}, Cchar, Ptr{Cuchar}, Ptr{Cint}))
    myevalhl_c = @cfunction($myevalhl, Cvoid, (Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cdouble, Ptr{Cdouble},
                                               Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Cint,
                                               Ptr{Cuchar}, Ptr{Cint}))
    myevalhlp_c = @cfunction($myevalhlp, Cvoid, (Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cdouble,
                                                 Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                                                 Ptr{Cuchar}, Ptr{Cint}))



    # additional problem description
    x = zeros(n)

    l = [-10.0 -1.0e20] # vector of lower bounds of x
    u = [10.0 1.0e20] # upper

    if equatn == nothing
        equatn = zeros(m)
    end
    equatn = map(Cuchar, equatn) # for each i ∈ {1, m}, equatn[i] == ? iff the mth constraint is an equation (instead of an inequality)
    linear = map(Cuchar, zeros(m)) # same, for linearity
    linear[1] = Cuchar(1)

    coded = map(Cuchar, zeros(11)) # same, for being effectively coded by the user
    if f != nothing
        coded[1] = Cuchar(1)
    end
    if g != nothing
        coded[2] = Cuchar(1)
    end
    if h != nothing
        coded[3] = Cuchar(1)
    end
    if c != nothing
        coded[4] = Cuchar(1)
    end
    if jac != nothing
        coded[5] = Cuchar(1)
    end
    if hc != nothing
        coded[6] = Cuchar(1)
    end

    lambda = zeros(m)

    checkder = 0 # 1 to check derivatives with finite differences

    jcnnzmax = 4 # max number of non nulls in the sparse Jacobian of the constraints
    hnnzmax = 7 # same, for the hessian of the objective function

    epsfeas = 1e-8 # epsilon for feasibility
    epsopt = 1e-8 #  epsilon for optimality
    efstain = sqrt(epsfeas) # epsilons for stopping at unfeasible stationary points
    eostain = epsopt^1.5 # I still don't understand these vars perfectly
    efacc = sqrt(epsfeas) # at which level of feasibility newton-like acceleration starts
    eoacc = sqrt(epsopt)

    outputfnm = "algencan.out"
    specfnm = ""
    nvparam = 1
    vparam = ["ITERATIONS-OUTPUT-DETAIL 10"]

    fval = cnorm = snorm = nlpsupn = inform = 0.0

    ccall(algencan_optimize, Cvoid, (
                            Ptr{Nothing},    # *myevalf,
                            Ptr{Nothing},    # *myevalg,
                            Ptr{Nothing},    # *myevalh,
                            Ptr{Nothing},    # *myevalc,
                            Ptr{Nothing},    # *myevaljac,
                            Ptr{Nothing},    # *myevalhc,
                            Ptr{Nothing},    # *myevalfc,
                            Ptr{Nothing},    # *myevalgjac,
                            Ptr{Nothing},    # *myevalgjacp,
                            Ptr{Nothing},    # *myevalhl,
                            Ptr{Nothing},    # *myevalhlp,
                            Cint,            # jcnnzmax,
                            Cint,            # hnnzmax,
                            Ref{Cdouble},    # *epsfeas,
                            Ref{Cdouble},    # *epsopt,
                            Ref{Cdouble},    # *efstain,
                            Ref{Cdouble},    # *eostain,
                            Ref{Cdouble},    # *efacc,
                            Ref{Cdouble},    # *eoacc,
                            Cstring,         # *outputfnm,
                            Cstring,         # *specfnm,
                            Cint,            # nvparam,
                            Ptr{Ptr{Cuchar}},# **vparam,
                            Cint,            # n,
                            Ref{Cdouble},    # *x,
                            Ref{Cdouble},    # *l,
                            Ref{Cdouble},    # *u,
                            Cint,            # m,
                            Ref{Cdouble},    # *lambda,
                            Ref{Cuchar},     # *equatn,
                            Ref{Cuchar},     # *linear,
                            Ref{Cuchar},     # *coded,
                            Cuchar,          # checkder,
                            Ref{Cdouble},    # *f,
                            Ref{Cdouble},    # *cnorm,
                            Ref{Cdouble},    # *snorm,
                            Ref{Cdouble},    # *nlpsupn,
                            Ref{Cint}        # *inform
                           ),
          myevalf_c, myevalg_c, myevalh_c, myevalc_c, myevaljac_c, myevalhc_c, myevalfc_c,
          myevalgjac_c, myevalgjacp_c, myevalhl_c, myevalhlp_c, jcnnzmax, hnnzmax,
          epsfeas, epsopt, efstain, eostain, efacc, eoacc, outputfnm, specfnm,
          nvparam, vparam, n, x, l, u, m, lambda,
          equatn, linear, coded, checkder, fval, cnorm, snorm,
          nlpsupn, inform
         )

    # OUTPUTS ARE NOT WORKING CORRECTLY (JULIA'S PROBABLY PASSING BY VALUE)

    return x, fval
end

function unload()
    Libdl.dlclose(algencan)
end

end
