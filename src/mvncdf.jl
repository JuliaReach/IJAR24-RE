###
#   This file is part of the ProbabilisticReachability package.
#
#   Definition of multivariate gaussian cdf for use in gaussian copula is described in:
#
#       https://discourse.julialang.org/t/mvn-cdf-have-it-coded-need-help-getting-integrating-into-distributions-jl/38631/16
#
#   See accompanying license in the `qsimvnv` function below.
###

function qsimvnv(Σ, a, b; m=nothing)
    #= rev 1.13
    	This function uses an algorithm given in the paper
    	"Numerical Computation of Multivariate Normal Probabilities", in
    	J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
    	Alan Genz, WSU Math, PO Box 643113, Pullman, WA 99164-3113
    	Email : alangenz@wsu.edu
    	The primary references for the numerical integration are
    	"On a Number-Theoretical Integration Method"
    	H. Niederreiter, Aequationes Mathematicae, 8(1972), pp. 304-11, and
    	"Randomization of Number Theoretic Methods for Multiple Integration"
    	R. Cranley and T.N.L. Patterson, SIAM Journal on Numerical Analysis, 13(1976), pp. 904-14.

    	Re-coded in Julia from the MATLAB function qsimvnv(m,r,a,b)

    	Alan Genz is the author the MATLAB qsimvnv() function.
    	Alan Genz software website: http://archive.is/jdeRh
    	Source code to MATLAB qsimvnv() function: http://archive.is/h5L37
    	% QSIMVNV(m,r,a,b) and _chlrdr(r,a,b)
    	%
    	% Copyright (C) 2013, Alan Genz,  All rights reserved.
    	%
    	% Redistribution and use in source and binary forms, with or without
    	% modification, are permitted provided the following conditions are met:
    	%   1. Redistributions of source code must retain the above copyright
    	%      notice, this list of conditions and the following disclaimer.
    	%   2. Redistributions in binary form must reproduce the above copyright
    	%      notice, this list of conditions and the following disclaimer in
    	%      the documentation and/or other materials provided with the
    	%      distribution.
    	%   3. The contributor name(s) may not be used to endorse or promote
    	%      products derived from this software without specific prior
    	%      written permission.
    	% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    	% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    	% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
    	% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
    	% COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
    	% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
    	% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
    	% OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    	% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
    	% TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF USE
    	% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    	%

    	Julia dependencies
    	Distributions
    	PDMats
    	Primes
    	Random
    	LinearAlgebra

    	=#

    if isnothing(m)
        m = 1000 * size(Σ, 1)  # default is 1000 * dimension
    end

    # check for proper dimensions
    n = size(Σ, 1)
    nc = size(Σ, 2) # assume square Cov matrix nxn
    # check dimension > 1
    n >= 2 || throw(ErrorException("dimension of Σ must be 2 or greater. Σ dimension: $(size(Σ))"))
    n == nc || throw(DimensionMismatch("Σ matrix must be square. Σ dimension: $(size(Σ))"))

    # check dimensions of lower vector, upper vector, and cov matrix match
    (n == size(a, 1) == size(b, 1)) ||
        throw(DimensionMismatch("iconsistent argument dimensions. Sizes: Σ $(size(Σ))  a $(size(a))  b $(size(b))"))

    # check that a and b are column vectors; if row vectors, fix it
    if size(a, 1) < size(a, 2)
        a = transpose(a)
    end
    if size(b, 1) < size(b, 2)
        b = transpose(b)
    end

    # check that lower integration limit a < upper integration limit b for all elements
    all(a .<= b) ||
        throw(ArgumentError("lower integration limit a must be <= upper integration limit b"))

    # check that Σ is positive definite; if not, print warning
    isposdef(Σ) || @warn "covariance matrix Σ fails positive definite check"

    # check if Σ, a, or b contains NaNs
    if any(isnan.(Σ)) || any(isnan.(a)) || any(isnan.(b))
        p = NaN
        e = NaN
        return (p, e)
    end

    # check if a==b
    if a == b
        p = 0.0
        e = 0.0
        return (p, e)
    end

    # check if a = -Inf & b = +Inf
    if all(a .== -Inf) && all(b .== Inf)
        p = 1.0
        e = 0.0
        return (p, e)
    end

    # check input Σ, a, b are floats; otherwise, convert them
    if eltype(Σ) <: Signed
        Σ = float(Σ)
    end

    if eltype(a) <: Signed
        a = float(a)
    end

    if eltype(b) <: Signed
        b = float(b)
    end

    ##################################################################
    #
    # Special cases: positive Orthant probabilities for 2- and
    # 3-dimesional Σ have exact solutions. Integration range [0,∞]
    #
    ##################################################################

    if all(a .== zero(eltype(a))) && all(b .== Inf) && n <= 3
        Σstd = sqrt.(diag(Σ))
        Rcorr = cov2cor(Σ, Σstd)

        if n == 2
            p = 1 / 4 + asin(Rcorr[1, 2]) / (2π)
            e = eps()
        elseif n == 3
            p = 1 / 8 + (asin(Rcorr[1, 2]) + asin(Rcorr[2, 3]) + asin(Rcorr[1, 3])) / (4π)
            e = eps()
        end

        return (p, e)
    end

    ##################################################################
    #
    # get lower cholesky matrix and (potentially) re-ordered integration vectors
    #
    ##################################################################

    (ch, as, bs) = _chlrdr(Σ, a, b) # ch =lower cholesky; as=lower vec; bs=upper vec

    ##################################################################
    #
    # quasi-Monte Carlo integration of MVN integral
    #
    ##################################################################

    ### setup initial values
    ai = as[1]
    bi = bs[1]
    ct = ch[1, 1]

    unitnorm = Normal() # unit normal distribution
    rng = RandomDevice()

    # if ai is -infinity, explicitly set c=0
    # implicitly, the algorithm classifies anything > 9 std. deviations as infinity
    if ai > -9 * ct
        if ai < 9 * ct
            c1 = Distributions.cdf.(unitnorm, ai / ct)
        else
            c1 = 1.0
        end
    else
        c1 = 0.0
    end

    # if bi is +infinity, explicitly set d=0
    if bi > -9 * ct
        if bi < 9 * ct
            d1 = Distributions.cdf(unitnorm, bi / ct)
        else
            d1 = 1.0
        end
    else
        d1 = 0.0
    end

    #n=size(Σ,1) 	# assume square Cov matrix nxn
    cxi = c1# initial cxi; genz uses ci but it conflicts with Lin. Alg. ci variable
    dci = d1 - cxi# initial dcxi
    p = 0.0# probability = 0
    e = 0.0# error = 0

    # Richtmyer generators
    ps = sqrt.(primes(Int(floor(5 * n * log(n + 1) / 4)))) # Richtmyer generators
    q = ps[1:(n - 1), 1]
    ns = 12
    nv = Int(max(floor(m / ns), 1))

    Jnv    = ones(1, nv)
    cfill  = transpose(fill(cxi, nv)) # evaluate at nv quasirandom points row vec
    dpfill = transpose(fill(dci, nv))
    y      = zeros(n - 1, nv)# n-1 rows, nv columns, preset to zero

    #=Randomization loop for ns samples
    	j is the number of samples to integrate over,
    		but each with a vector nv in length
    	i is the number of dimensions, or integrals to compute =#

    for j in 1:ns# loop for ns samples
        c  = copy(cfill)
        dc = copy(dpfill)
        pv = copy(dpfill)
        for i in 2:n
            x = transpose(abs.(2.0 .* mod.((1:nv) .* q[i - 1] .+ rand(rng), 1) .- 1)) # periodizing transformation
            # note: the rand() is not broadcast -- it's a single random uniform value added to all elements
            y[i - 1, :] = quantile.(unitnorm, c .+ x .* dc)
            s = transpose(ch[i, 1:(i - 1)]) * y[1:(i - 1), :]
            ct = ch[i, i]# ch is cholesky matrix
            ai = as[i] .- s
            bi = bs[i] .- s
            c = copy(Jnv)# preset to 1 (>9 sd, +∞)
            d = copy(Jnv)# preset to 1 (>9 sd, +∞)

            c[findall(x -> isless(x, -9 * ct), ai)] .= 0.0# If < -9 sd (-∞), set to zero
            d[findall(x -> isless(x, -9 * ct), bi)] .= 0.0# if < -9 sd (-∞), set to zero
            tstl = findall(x -> isless(abs(x), 9 * ct), ai)# find cols between -9 and +9 sd (-∞ to +∞)
            c[tstl] .= Distributions.cdf.(unitnorm, ai[tstl] / ct)# for those, compute Normal CDF
            tstl = (findall(x -> isless(abs(x), 9 * ct), bi))# find cols between -9 and +9 sd (-∞ to +∞)
            d[tstl] .= Distributions.cdf.(unitnorm, bi[tstl] / ct)
            @. dc = d - c
            @. pv = pv * dc
        end # for i=
        d = (mean(pv) - p) / j
        p += d
        e = (j - 2) * e / j + d^2
    end # for j=

    e = 3 * sqrt(e) # error estimate is 3 times standard error with ns samples

    return (p, e)  # return probability value and error estimate
end # function qsimvnv

function _chlrdr(Σ, a, b)
    # Rev 1.13

    # define constants
    # 64 bit machine error 1.0842021724855e-19 ???
    # 32 bit machine error 2.220446049250313e-16 ???
    ep = 1e-10 # singularity tolerance
    if Sys.WORD_SIZE == 64
        fpsize = Float64
        ϵ = eps(0.0) # 64-bit machine error
    else
        fpsize = Float32
        ϵ = eps(0.0f0) # 32-bit machine error
    end

    if !@isdefined sqrt2π
        sqrt2π = √(2π)
    end

    # unit normal distribution
    unitnorm = Normal()

    n = size(Σ, 1) # covariance matrix n x n square

    ckk = 0.0
    dem = 0.0
    am = 0.0
    bm = 0.0
    ik = 0.0

    if eltype(Σ) <: Signed
        c = copy(float(Σ))
    else
        c = copy(Σ)
    end

    if eltype(a) <: Signed
        ap = copy(float(a))
    else
        ap = copy(a)
    end

    if eltype(b) <: Signed
        bp = copy(float(b))
    else
        bp = copy(b)
    end

    d = sqrt.(diag(c))
    for i in 1:n
        if d[i] > 0.0
            c[:, i] /= d[i]
            c[i, :] /= d[i]
            ap[i] = ap[i] / d[i]     # ap n x 1 vector
            bp[i] = bp[i] / d[i]     # bp n x 1 vector
        end
    end

    y = zeros(fpsize, n) # n x 1 zero vector to start

    for k in 1:n
        ik = k
        ckk = 0.0
        dem = 1.0
        s = 0.0
        #pprinta(c)
        for i in k:n
            if c[i, i] > ϵ  # machine error
                cii = sqrt(max(c[i, i], 0))

                if i > 1 && k > 1
                    s = (c[i, 1:(k - 1)] .* y[1:(k - 1)])[1]
                end

                ai = (ap[i] - s) / cii
                bi = (bp[i] - s) / cii
                de = Distributions.cdf(unitnorm, bi) - Distributions.cdf(unitnorm, ai)

                if de <= dem
                    ckk = cii
                    dem = de
                    am = ai
                    bm = bi
                    ik = i
                end
            end # if c[i,i]> ϵ
        end # for i=
        i = n

        if ik > k
            ap[ik], ap[k] = ap[k], ap[ik]
            bp[ik], bp[k] = bp[k], bp[ik]

            c[ik, ik] = c[k, k]

            if k > 1
                c[ik, 1:(k - 1)], c[k, 1:(k - 1)] = c[k, 1:(k - 1)], c[ik, 1:(k - 1)]
            end

            if ik < n
                c[(ik + 1):n, ik], c[(ik + 1):n, k] = c[(ik + 1):n, k], c[(ik + 1):n, ik]
            end

            if k <= (n - 1) && ik <= n
                c[(k + 1):(ik - 1), k], c[ik, (k + 1):(ik - 1)] = transpose(c[ik, (k + 1):(ik - 1)]),
                                                                  transpose(c[(k + 1):(ik - 1), k])
            end
        end # if ik>k

        if ckk > k * ep
            c[k, k] = ckk
            if k < n
                c[k:k, (k + 1):n] .= 0.0
            end

            for i in (k + 1):n
                c[i, k] /= ckk
                c[i:i, (k + 1):i] -= c[i, k] * transpose(c[(k + 1):i, k])
            end

            if abs(dem) > ep
                y[k] = (exp(-am^2 / 2) - exp(-bm^2 / 2)) / (sqrt2π * dem)
            else
                if am < -10
                    y[k] = bm
                elseif bm > 10
                    y[k] = am
                else
                    y[k] = (am + bm) / 2
                end
            end # if abs
        else
            c[k:n, k] .== 0.0
            y[k] = 0.0
        end # if ckk>ep*k
    end # for k=

    return (c, ap, bp)
end # function _chlrdr

#####################
# CDF for reach-sets
#####################

import ProbabilityBoundsAnalysis: cdf
const IA = IntervalArithmetic

_check_time(R, t::Number) = t ∈ tspan(R)
_check_time(R, t::IntervalArithmetic.Interval) = t ⊆ tspan(R)

# focal_elements: vector of interval boxes or intervals
# masses: vector of reals
function get_pbox(R::TaylorModelReachSet{N}, t::Number, focal_elems::AbstractVector{IA.Interval{N}},
                  masses::Vector{Float64}=ones(length(focal_elems)); check::Bool=true) where {N}
    nelems = length(focal_elems)

    if check
        _check_time(R, t) || return pbox()
    end

    out = Vector{IA.Interval{N}}(undef, nelems)
    @inbounds for i in 1:nelems
        Bi = IntervalBox(focal_elems[i])
        out[i] = evaluate(R, t, Bi)[1]
    end

    return pbox(out, masses)
end

function get_pbox(R::TaylorModelReachSet{N}, t::Number,
                  focal_elems::Vector{IA.IntervalBox{T,Float64}},
                  masses::Vector{Float64}=ones(length(focal_elems)); check::Bool=true) where {N,T}
    nelems = length(focal_elems)

    if check
        _check_time(R, t) || return pbox()
    end

    out = Vector{IA.IntervalBox{LazySets.dim(R),Float64}}(undef, nelems)
    @inbounds for i in 1:nelems
        Bi = focal_elems[i]
        out[i] = evaluate(R, t, Bi)
    end

    outPboxes = Vector{pbox}(undef, LazySets.dim(R))

    for i in 1:LazySets.dim(R)
        theseIntervals = [out[j][i] for j in 1:nelems]
        outPboxes[i] = pbox(theseIntervals, masses)
    end

    return outPboxes
end

function get_pbox(R::AbstractVector{<:TaylorModelReachSet{N}}, t::Number,
                  focal_elems::AbstractVector{IA.Interval{N}},
                  masses::Vector{Float64}=ones(length(focal_elems))) where {N}
    return env((get_pbox(Ri, t, focal_elems, masses) for Ri in R)...)
end

function get_pbox(R::AbstractVector{<:TaylorModelReachSet{N}}, dt::IntervalArithmetic.Interval,
                  focal_elems::AbstractVector{IA.Interval{N}},
                  masses::Vector{Float64}=ones(length(focal_elems))) where {N}
    idx = findall(Ri -> !isempty(tspan(Ri) ∩ dt), R)
    return env((get_pbox(Ri, tspan(Ri) ∩ dt, focal_elems, masses; check=false) for Ri in
                                                                                   view(R, idx))...)
end

function get_pbox(R::AbstractVector{<:TaylorModelReachSet{N}}, t, pb::pbox) where {N}
    focal_elems, masses = split(pb)
    return get_pbox(R, t, focal_elems, masses)
end

function get_pbox(R::TaylorModelReachSet{N}, t, pb::pbox) where {N}
    focal_elems, masses = split(pb)
    return get_pbox(R, t, focal_elems, masses)
end

function get_pbox(R::TaylorModelReachSet{N}, t, pb::nPbox; quick=false) where {N}
    if quick
        return get_pbox_quick(R, t, pb)
    end
    focal_elems, masses = split(pb)
    return get_pbox(R, t, focal_elems, masses)
end

function get_pbox_quick(R::TaylorModelReachSet{N}, t, pb::nPbox) where {N}
    pb = ProbabilisticReachability.normalize(pb)
    margs = pb.marginals

    ranges = [range(mar) for mar in margs]

    outPboxes = Vector{pbox}(undef, LazySets.dim(R))

    for i in 1:LazySets.dim(R)
        fes, masses = split(margs[i])

        if i == 1
            Bi = fes .× Ref(IntervalBox(ranges[2:end]))
        elseif i == LazySets.dim(R)
            Bi = Ref(IntervalBox(ranges[1:(end - 1)])) .× fes
        else
            lefts = IntervalBox(ranges[1:(i - 1)])
            rights = IntervalBox(ranges[(i + 1):end])
            Bi = Ref(lefts) .× fes .× Ref(rights)
        end
        out = evaluate.(Ref(R), t, Bi)
        theseIntervals = [out[j][i] for j in 1:length(out)]
        outPboxes[i] = pbox(theseIntervals, masses)
    end

    return outPboxes
end

function get_pbox(R::AbstractVector{<:TaylorModelReachSet{N}}, t, pb::nPbox; quick=false) where {N}
    focal_elems, masses = split(pb)
    return get_pbox(R, t, focal_elems, masses; quick=quick)
end

function get_pbox(sol::ReachabilityAnalysis.AbstractSolution, t,
                  focal_elems::AbstractVector{IA.Interval{N}},
                  masses::Vector{Float64}=ones(length(focal_elems))) where {N}
    return get_pbox(sol(t), t, focal_elems, masses)
end

function get_pbox(sol::ReachabilityAnalysis.AbstractSolution, t, pb::pbox; check_range::Bool=true)
    if check_range
        pb0 = get_initial_pbox(sol)
        pb0r = range(pb0)
        pbr = range(pb)
        pb0r == pbr ||
            throw(ArgumentError("the range of the given p-box, $pbr, does not correspond to the range of the solutions's p-box, $pb0r"))
    end
    pb_norm = normalize(pb)
    return get_pbox(sol(t), t, pb_norm)
end

function get_pbox(sol::ReachabilityAnalysis.AbstractSolution, t, pb::nPbox; check_range::Bool=true,
                  quick=false)
    if check_range
        pb0 = get_initial_pbox(sol)
        pb0r = range(pb0)
        pbr = range(pb)
        pb0r == pbr ||
            throw(ArgumentError("the range of the given p-box, $pbr, does not correspond to the range of the solutions's p-box, $pb0r"))
    end
    pb_norm = normalize(pb)
    return get_pbox(sol(t), t, pb_norm; quick=quick)
end

function get_pbox(sol::ReachabilityAnalysis.AbstractSolution, t; quick=false)
    pb = get_initial_pbox(sol; normalized=true)
    if typeof(pb) == pbox
        return get_pbox(sol(t), t, pb)
    end
    return get_pbox(sol(t), t, pb; quick=quick)
end
