######
# This file is part of the ProbabilisticReachability package
#
#   Definition of N-D copula and constructors
#
######

##
#   A multivariate copula
##
struct nCopula{N}

    #cdf :: Union{Array{Float64, N}, Missing}
    #Func :: Union{<:Function, Missing}

    Func::Function
    fam::String
    param::Any
end

# Return copula cdf for a vector u
function (obj::nCopula)(u::Vector{<:Real})
    if ismissing(obj.Func)
        return StepInterp(obj, u)
    end
    return obj.Func(u)
end

cdf(C::nCopula, u::Vector{<:Real}) = C(u)

dims(X::nCopula{N}) where {N} = N

###
#   Returns interval value for copula which has be evaluated on a grid.
#   Called when the function of C is not known. Will not yet be considered
###
function StepInterp(C::nCopula, u::Array{Float64,1})
    throw(ArgumentError("Step interpolation not yet implemented"))
end

#   Returns gaussian copula cdf at u, for correlation matrix Σ
function GauCopCDF(Σ, u)
    xs = quantile.(Normal(), u)
    a = ones(length(xs)) * -Inf

    return qsimvnv(Σ, a, xs)[1]
end

#   π/indepdenent/product copula function
indep(u) = prod(u)

###
#   Parametric copula constructors
###
IndepCopula(Ndim::Int=2) = nCopula{Ndim}(indep, "π", missing)

function GaussianCopula(Σ)
    posDef = isposdef(Σ)
    if !posDef
        throw(ArgumentError("Not valid correlation matrix"))
    end

    dims = size(Σ)[1]

    GauFunc(u) = GauCopCDF(Σ, u)
    return nCopula{dims}(GauFunc, "Gaussian", Σ)
end

###
#   Perfect and opposite copulas
###

function perfect(Ndim::Int=2)
    Func(u) = minimum(u)
    return nCopula{Ndim}(Func, "M", missing)
end

# opposite copula only defined for 2D
function opposite()
    Func(u) = max(u[1] + u[2] - 1, 0)
    return nCopula{2}(Func, "W", missing)
end

###
#   Sampling function for
###

function rand(C::nCopula, nsamples::Int)
    if C.fam == "Gaussian"
        return rand_gau_copula(C.param, nsamples)
    end

    if C.fam == "π"
        return rand(nsamples, dims(C))
    end

    throw(ArgumentError("Copula sampler not implemented for family $(C.fam)"))
end

function rand_gau_copula(Σ, nsamples)
    L = cholesky(Σ).L

    Z = rand(Normal(), nsamples, size(Σ, 2))
    X = transpose(L * transpose(Z))

    return cdf.(Normal(), X)
end
