module IntervalStatistic

include("Average.jl")
include("Variance.jl")

function mean(values :: Vector, method :: Average.BaseEstimator)
    Average.estimate(values, method)
end

function var(values :: Vector, method :: Variance.BaseEstimator)
    Variance.estimate(values, method)
end

end
