module Check
using StatsBase
using Distributions
using ValidatedNumerics

abstract BaseCheck
function isDistribution(values :: Vector, check::BaseCheck)
    error("Unsupported check $check")
end


abstract NormalCheck <: BaseCheck

immutable ChiSquareCheck <: BaseCheck
    error_probability :: Real
    distribution
    calculateCoefficient :: Function
end

function SimpleChiSquareCheck(error_probability, distribution)
    ChiSquareCheck(error_probability, distribution, (n) -> round(Int, n / 5))
end

function isDistribution(values :: Vector, check :: ChiSquareCheck)
    values_count = length(values)
    intervals_count = check.calculateCoefficient(values_count)
    hist = fit(Histogram, values, nbins=intervals_count)
    n_probability = values_count/intervals_count

    chi_square = mapreduce(
        (x) -> chiSquareEstimate(x, values_count, check), +,
        zip(edgesToIntervals(hist.edges[1]), hist.weights)
    )

    freedom_degree_number = intervals_count - 1
    chi_square_quantile = getChiSquareQuantile(
        check.error_probability, freedom_degree_number
    )
    chi_square <= chi_square_quantile
end

function getChiSquareQuantile(value, freedom_degree_number)
    d = Chisq(freedom_degree_number)
    quantile(d, value)
end

function edgesToIntervals(edges)
    collected = collect(edges)
    downs = collected[1:end-1]
    ups = collected[2:end]
    [@interval(i[1], i[2]) for i in zip(downs, ups)]
end

function chiSquareEstimate(edgeAndValue, total_count, check :: ChiSquareCheck)
    interval = edgeAndValue[1]
    count_in_interval = edgeAndValue[2]
    d = check.distribution
    theoretical_probability = cdf(d, interval.hi) - cdf(d, interval.lo)
    theoretical_hist = total_count*theoretical_probability
    ((count_in_interval - theoretical_hist)^2)/theoretical_hist
end

end
