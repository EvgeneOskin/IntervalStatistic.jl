module Check
using StatsBase
using Distributions
using Debug

abstract BaseCheck
function isDistribution(values :: Vector, check::BaseCheck)
    error("Unsupported check $check")
end


abstract NormalCheck <: BaseCheck

immutable ChiSquareCheck <: BaseCheck
    error_probability :: Real
    calculateCoefficient :: Function
end
function SimpleChiSquareCheck(error_probability)
    ChiSquareCheck(error_probability, (n)->5)
end

@debug function isDistribution(values :: Vector, check :: ChiSquareCheck)
    values_count = size(values, 1)
    intervals = generateIntervals(values, check)
    hist = fit(Histogram, values, intervals)
    intervals_count = size(intervals, 1)
    n_probability = values_count/intervals_count

    chi_square = n_probability*mapreduce(
        (x) -> (x/n_probability - 1)^2, +, hist.weights
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

function generateIntervals(values :: Vector, method :: ChiSquareCheck)
    count = size(values, 1)
    interval_count = count/method.calculateCoefficient(count)
    max_values, min_value = maximum(values), minimum(values)
    step = (max_values - min_value) / interval_count
    min_value:step:(max_values + step*0.01)
end

end
