module Check
using StatsBase
using Distributions

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
    ChiSquareCheck(error_probability, (n) -> Int(round(n / 5)))
end

function isDistribution(values :: Vector, check :: ChiSquareCheck)
    values_count = size(values, 1)
    intervals_count = check.calculateCoefficient(values_count)
    hist = fit(Histogram, values, nbins=intervals_count)
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

end
