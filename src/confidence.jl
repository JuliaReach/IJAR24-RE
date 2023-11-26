# ======================
# Confidence intervals
# ======================

# return the p-value confidence interval (p between 0 an 1)
# for the given p-box x
function confidence(x::pbox, p)
    p_lo = (1 - p) / 2
    p_hi = 1 - p_lo

    c_lo = cut(x, p_lo)
    c_hi = cut(x, p_hi)

    return interval(c_lo.lo, c_hi.hi)
end

function confidence(x::nPbox, p)
    focal_elems, poss = make_cons(x)

    cs = 1 .- poss
    D = ProbabilisticReachability.dims(x)
    fe_sets = Vector{IA.IntervalBox{D,Float64}}(undef, length(p))
    for (i, ip) in enumerate(p)
        id = findfirst(cs .<= ip)
        fe_sets[i] = focal_elems[id]
    end
    return fe_sets
end

function confidence(sol::ReachSolution, p; pbinit=get_initial_pbox(sol),
                    normalize=true, cons=true)
    if LazySets.dim(sol) == 1
        if cons
            confidence_univariate_cons(sol, p, pbinit; normalize=normalize)
        else
            confidence_univariate(sol, p, pbinit; normalize=normalize)
        end
    else
        if cons
            confidence_multivariate_cons(sol, p, pbinit; normalize=normalize)
        else
            confidence_multivariate(sol, p, pbinit; normalize=normalize)
        end
    end
end

function confidence_univariate(sol::ReachSolution, p, pbinit; normalize=true)
    # extract the confidence intervals
    numprob = length(p)
    numr = length(sol)
    N = Float64
    IT = LazySets.Interval{N,IA.Interval{N}}
    RT = ReachSet{N,IT}
    FT = Flowpipe{N,RT,Vector{RT}}
    out = [Flowpipe(undef, RT, numr) for _ in 1:numprob]

    if normalize
        pbinit = ProbabilisticReachability.normalize(pbinit)
    end
    focal_elems, masses = split(pbinit)

    for (i, R) in enumerate(sol) # loop over the reach-set

        # TEMP should use tspan(R)
        tt = narrow(interval(tstart(R), tend(R)))

        pb = get_pbox(R, tt, focal_elems, masses; check=true)

        for (j, pj) in enumerate(p) # loop over the probabilities
            y = confidence(pb, pj)
            Y = LazySets.Interval(y)
            out[j].Xk[i] = ReachSet(Y, tt)
        end
    end

    return out
end

function confidence_multivariate(sol::ReachSolution, p, pbinit; normalize=true)
    # extract the confidence intervals
    numprob = length(p)
    numr = length(sol)
    N = Float64
    ST = Hyperrectangle{N,Vector{N},Vector{N}}
    RT = ReachSet{N,ST}
    FT = Flowpipe{N,RT,Vector{RT}}
    out = [Flowpipe(undef, RT, numr) for _ in 1:numprob]

    if normalize
        pbinit = ProbabilisticReachability.normalize(pbinit)
    end
    focal_elems, masses = split(pbinit)

    for (i, R) in enumerate(sol) # loop over the reach-set
        tt = narrow(interval(tstart(R), tend(R)))
        pb = get_pbox(R, tt, focal_elems, masses; check=true)

        for (j, pj) in enumerate(p) # loop over the probabilities
            y = IntervalBox([confidence(x, pj) for x in pb])
            Y = convert(Hyperrectangle, y)
            out[j].Xk[i] = ReachSet(Y, tt)
        end
    end

    return out
end

function confidence_univariate_cons(sol::ReachSolution, p, pbinit; normalize=true,
                                    partition=ones(Int, LazySets.dim(sol)),
                                    ntdiv=1)

    # extract the confidence intervals
    numprob = length(p)
    N = Float64
    IT = LazySets.Interval{N,IA.Interval{N}}
    RT = ReachSet{N,IT}
    FT = Flowpipe{N,RT,Vector{RT}}
    out = [Flowpipe(Vector{RT}()) for _ in 1:numprob]

    if normalize
        pbinit = ProbabilisticReachability.normalize(pbinit)
    end

    fe_sets = confidence.(pbinit, p)

    # each fe_sets should be contained in -1 .. 1
    fe_sets[isone.(p)] .= -1 .. 1

    for (j, pj) in enumerate(p) # loop over the probabilities
        Bi = IntervalBox(fe_sets[j])
        for (i, R) in enumerate(sol) # loop over the reach-set
            # first idea
            #tt = narrow(interval(tstart(R), tend(R)))
            #y =  evaluate(R, tt, Bi)
            #Y = LazySets.Interval(y[1])
            #out[j].Xk[i] = ReachSet(Y, tt)

            # second idea: union of boxes
            # list of ReachSet whose set rep is IntervalBox
            #Y = overapproximate(R, Hyperrectangle, partition=partition)
            #append!(out[j], Y)

            # second idea: single zonotope at given subdomain
            Y = overapproximate(R, Zonotope; dom=Bi)
            Y = ReachabilityAnalysis.reconstruct(Y, convert(LazySets.Interval, set(Y)))
            push!(out[j], Y)
        end
    end
    return out
end

function confidence_multivariate_cons(sol::ReachSolution, p, pbinit; normalize=true,
                                      partition=ones(Int, dim(sol)),
                                      ntdiv=1)
    # extract the confidence intervals
    numprob = length(p)
    N = Float64
    ST = Zonotope{N,Vector{N},Matrix{N}}
    RT = ReachSet{N,ST}
    FT = Flowpipe{N,RT,Vector{RT}}
    out = [Flowpipe(Vector{RT}()) for _ in 1:numprob]

    if normalize
        pbinit = ProbabilisticReachability.normalize(pbinit)
    end
    fe_sets = confidence(pbinit, p)

    fe_sets[isone.(p)] .= Ref(IntervalBox(-1 .. 1, dims(pbinit)))

    for (j, pj) in enumerate(p) # loop over the probabilities
        for (i, R) in enumerate(sol) # loop over the reach-set
            push!(out[j], overapproximate(R, Zonotope; dom=fe_sets[j]))
        end
    end
    return out
end

function mean_flowpipe(sol::ReachSolution)
    # extract the confidence intervals
    numr = length(sol)
    N = Float64
    IT = LazySets.Interval{N,IA.Interval{N}}
    RT = ReachSet{N,IT}
    FT = Flowpipe{N,RT,Vector{RT}}
    out = Flowpipe(undef, RT, numr)

    for (i, s) in enumerate(sol) # loop over the reach-set

        # temp, should use tspan(s)
        tt = interval(tstart(s) + 1e-14, tend(s) - 1e-14)

        pb = get_pbox(sol, tt)

        y = mean(pb)
        R = ReachSet(LazySets.Interval(y), tt)
        out.Xk[i] = R
    end

    return out
end
