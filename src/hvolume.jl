######
# This file is part of the ProbabilisticReachability package
#
#   H-volume & related functions. H-volume finds the probability measure in some box from a
#   multivariate cdf.
#
######

function Hvolume(measure::Union{nCopula,nPbox}, dom::IntervalBox)
    Ndims = length(dom)
    Nverts = 2^Ndims

    Js = [lefts(dom) rights(dom)]                  # Matrix of endpoint

    verts = vertices(dom)                           # All vertices of box

    if size(verts, 1) != Nverts
        throw(ArgumentError("Number of returned vertices does not match needed number"))
    end

    measureVals = zeros(Nverts, 1)

    for i in 1:Nverts
        measureVals[i] = measure(verts[i, :])    # Eval vertex in measure

        Ns = 0
        for j in 1:Ndims
            this = verts[i, j] == Js[j, 1]        # Find how many vertices are lower bounds
            Ns = Ns + this
        end

        sign = 1
        if mod(Ns, 2) == 1                           # If number of lower bounds are odd
            sign = -1                              # measure value will be subtracted
        end
        measureVals[i] = sign * measureVals[i]
    end

    return sum(measureVals)
end

##
#   Returns all iterations of endpoints from an interval box
##
function vertices(box::IntervalBox)
    Ndims = length(box)
    Nverts = 2^Ndims

    verts = zeros(Nverts, Ndims)

    for i in 1:Ndims
        A = [box[i].lo box[i].hi]
        j1 = 2^(i - 1)
        j2 = 2^(Ndims - i)

        As = repeat(A; outer=(j1, j2))

        verts[:, i] = As[:]
    end

    return verts
end

function HvolumeFast(cdfMatrix::Array{Float64,N}) where {N}
    sizes = size(cdfMatrix)

    Ndims = N
    Nverts = 2^Ndims

    dom = IntervalBox(ones(Ndims) .* [interval(0, 1)])
    Js = [lefts(dom) rights(dom)]                  # Matrix of endpoint

    verts = vertices(dom)                           # All vertices of box

    vals = Vector{typeof(cdfMatrix)}(undef, Nverts)

    R1s = [2:s for s in sizes]
    R2s = [1:(s - 1) for s in sizes]

    for i in 1:Nverts
        ThisR = deepcopy(R1s)
        for j in 1:Ndims
            if verts[i, j] == 0
                ThisR[j] = R2s[j]
            end
        end

        vals[i] = cdfMatrix[ThisR...]

        Ns = 0
        for j in 1:Ndims
            this = verts[i, j] == Js[j, 1]
            Ns = Ns + this
        end

        sign = 1
        if mod(Ns, 2) == 1
            sign = -1
        end
        vals[i] = sign * vals[i]
    end

    return sum(vals)
end

function make_cons(X::nPbox)
    margs = X.marginals
    C = X.copula

    n_dims = length(margs)

    ns = getfield.(margs, :n)
    alleq = all(ns .== ns[1])

    if !alleq
        throw(ArgumentError("All the steps of marginals must be equal"))
    end

    nn = ns[1]
    n_2 = Int(nn / 2)

    marg_cuts = Matrix{Interval{Float64}}(undef, n_2, n_dims)
    for j in 1:n_dims
        marg_cuts[:, j] = interval.(margs[j].u[1:n_2], reverse(margs[j].d[(n_2 + 1):end]))
    end

    joint_cuts = [IntervalBox(marg_cuts[j, :]) for j in 1:n_2]

    ps = range(0, 1; length=nn + 1)

    iis = ps[1:end]
    jjs = reverse(ps[2:end])

    pps = interval.(iis[1:n_2], jjs[1:n_2])

    joint_ps = IntervalBox.(pps, n_dims)

    probs = Hvolume.(Ref(C), joint_ps)

    poss = 1 .- probs

    ## Numerical instability of the mvgaussian integration gives some possibility values > 1

    index = findfirst(poss .>= 1)
    if !isnothing(index)
        poss[index:end] .= 1
        ln = length(poss)

        # removes all alpha cuts which are after poss >= 1
        [popat!(poss, length(poss)) for i in 1:(ln - index)]
        [popat!(joint_cuts, length(joint_cuts)) for i in 1:(ln - index)]
    end

    # Removes non-sorted alpha cuts
    if !issorted(poss)
        joint_cuts, poss = makesort(joint_cuts, poss)
    end

    return joint_cuts, poss
end

function make_cons(X::pbox)
    ns = X.n
    n_2 = Int(ns / 2)

    cuts = interval.(X.u[1:n_2], reverse(X.d[(n_2 + 1):end]))
    poss = range(0, 1; length=n_2 + 1)[1:(end - 1)]
    return cuts, collect(poss)
end

###
#   Converts a mv p-box or p-box in to a consonant belief function
###
function cons_split(X::PBOX)
    fe, poss = make_cons(X)

    is = [poss; 1]
    masses = is[2:end] .- is[1:(end - 1)]   # mass assignment is the difference between the possibility levels
    return fe, masses
end

###
#   Will remove alpha cuts which are not sorted. Retains riggor.
###
function makesort(cuts, poss)
    if issorted(poss)
        return cuts, poss
    end

    id = findfirst(poss[2:end] .< poss[1:(end - 1)])

    popat!(poss, id)
    popat!(cuts, id)
    return makesort(cuts, poss)
end
