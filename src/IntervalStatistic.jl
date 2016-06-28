module IntervalStatistic

include("Mean.jl")
include("Variance.jl")
include("Check.jl")

function mean(values :: Vector, method :: Mean.BaseEstimator)
    Mean.estimate(values, method)
end

function var(values :: Vector, method :: Variance.BaseEstimator)
    Variance.estimate(values, method)
end

function isDistribution(values :: Vector, method :: Check.BaseCheck)
    Check.isDistribution(values, method)
end

end
