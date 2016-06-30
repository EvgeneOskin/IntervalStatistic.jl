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

function LargeNChiSquareCheck(error_probability, distribution)
    # Mann H.B., Wald A., On the choice of number of intervals in the applciation of the chi-square test – 1942
    # Hollander M., Proshan F., Testing whether new is better than used – 1972
    ChiSquareCheck(error_probability, distribution, (n) -> calculateIntervalCountForLargeN(n))
end

function DahiyaChiSquareCheck(error_probability, distribution)
    # Dahiya R.C., Gurland J., How many classes in the Pearson chi-square test? – 1973
    ChiSquareCheck(error_probability, distribution, (n) -> round(Int, 1 + 3.32*log10(n)))
end


function calculateIntervalCountForLargeN(n)
    if n >= 500
        result = 4 * 0.75^0.2 * ( n -1 )^0.4
        round(Int, result)
    else
        error("$n show be greater then or equal 500.")
    end
end

function counts_and_histogram(values :: Vector, check :: ChiSquareCheck)
    values_count = length(values)
    intervals_count = check.calculateCoefficient(values_count)
    hist = fit(Histogram, values, nbins=intervals_count)
    intervals_and_weights = zip(edgesToIntervals(hist.edges[1]), hist.weights)
    (values_count, intervals_count, intervals_and_weights)
end


function histogram(values :: Vector, check :: ChiSquareCheck)
    _, _, hist = counts_and_histogram(values, check)
    hist
end

function isDistribution(values :: Vector, check :: ChiSquareCheck)
    values_count, intervals_count, hist = counts_and_histogram(values, check)
    n_probability = values_count/intervals_count

    chi_square = mapreduce(
        (x) -> chiSquareEstimate(x, values_count, check), +,
        hist
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
