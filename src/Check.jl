module Check
using StatsBase
using Distributions
using ValidatedNumerics

abstract BaseCheck

function isDistribution(values :: Vector, check::BaseCheck)
    error("Unsupported check $check")
end

abstract NormalCheck <: BaseCheck

abstract ChiSquareCheck <: BaseCheck

abstract ChiSquareCheckByIntervalCount <: ChiSquareCheck

abstract ChiSquareCheckByIntervalStep <: ChiSquareCheck

immutable LargeNChiSquareCheck <: ChiSquareCheckByIntervalCount
    # Mann H.B., Wald A., On the choice of number of intervals in the applciation of the chi-square test – 1942
    # Hollander M., Proshan F., Testing whether new is better than used – 1972
    error_probability :: Real
    distribution
end

function calculateBinCount(values :: Vector, check :: LargeNChiSquareCheck)
    n = length(values)
    if n >= 500
        result = 4 * 0.75^0.2 * ( n - 1 )^0.4
        round(Int, result)
    else
        error("$n show be greater then or equal 500.")
    end
end

immutable SturgesChiSquareCheck <: ChiSquareCheckByIntervalCount
    # Herbert A. Sturges, The Choice of a Class Interval – 1926
    error_probability :: Real
    distribution
end

function calculateBinCount(values :: Vector, check :: SturgesChiSquareCheck)
    n = length(values)
    round(Int, 1 + log2(n))
end

immutable ScottChiSquareCheck <: ChiSquareCheckByIntervalStep
    # David W. Scott, On optimal and data-based histograms –  1979
    error_probability :: Real
    distribution
end

function calculateBinStep(values :: Vector, check :: ScottChiSquareCheck)
    n = length(values)
    s = std(check.distribution)
    3.49 * s * n^(-1/3)
end

immutable TaylorChiSquareCheck <: ChiSquareCheckByIntervalStep
    # Charles C. Taylor, Akaike's information criterion and the histogram – 1987
    error_probability :: Real
    distribution
end

function calculateBinStep(values :: Vector, check :: TaylorChiSquareCheck)
    n = length(values)
    s = std(check.distribution)
    2.72 * s * n^(-1/3)
end

immutable FreedmanDiaconisChiSquareCheck <: ChiSquareCheckByIntervalStep
    # David Freedman, Persi Diaconis, On the Histogram as a Density Estimator: L2 Theory – 1981
    error_probability :: Real
    distribution
end

function calculateBinStep(values :: Vector, check :: FreedmanDiaconisChiSquareCheck)
    n = length(values)
    q1 = quantile(check.distribution, 0.25)
    q3 = quantile(check.distribution, 0.75)
    iqr = q3 - q1
    2 * iqr * n^(-1/3)
end

immutable DoaneChiSquareCheck <: ChiSquareCheckByIntervalCount
    # David P. Doane, Aesthetic Frequance Classifications – 2012
    error_probability :: Real
    distribution
end

function calculateBinCount(values :: Vector, check :: DoaneChiSquareCheck)
    average = mean(values)
    s = std(values)
    n = length(values)

    reduced_per_pow = (p) -> mapreduce(x -> (x - average)^p, +, values)
    sqrt_q = reduced_per_pow(3) / reduced_per_pow(2)^(3/2)
    sigma_sqrt_q = sqrt((6 * (n - 2)) / ((n + 1) * (n + 3)))

    result = 1 + log2(n) + log2(1 + sqrt_q/sigma_sqrt_q)
    round(Int, result)
end

immutable WichardChiSquareCheck <: ChiSquareCheckByIntervalCount
    # Jorg D. Wichard, Ronald Kuhne, Binding site detection via mutual information - 2008
    error_probability :: Real
    distribution
end

function calculateBinCount(values :: Vector, check :: WichardChiSquareCheck)
    average = mean(values)
    variance = var(values)
    n = length(values)

    reduced_per_pow = (p) -> mapreduce(x -> (x - average)^p, +, values)
    kurtosis = reduced_per_pow(4) / ((n - 1) * variance^2)
    result = 1 + log2(n) + log2(1 + kurtosis / sqrt(n / 6))
    round(Int, result)
end

function countsAndHistogram(values :: Vector, check :: ChiSquareCheck)
    values_count = length(values)
    intervals_count = calculateBinCount(values, check)
    hist = fit(Histogram, values, nbins=intervals_count)
    intervals_and_weights = zip(edgesToIntervals(hist.edges[1]), hist.weights)
    (values_count, intervals_count, intervals_and_weights)
end

function calculateBinCount(values :: Vector, check :: ChiSquareCheckByIntervalStep)
    values_count = length(values)
    max, min = maximum(values), minimum(values)
    intervals_step = calculateBinStep(values, check)
    round(Int, (max - min) / intervals_step)
end

function histogram(values :: Vector, check :: ChiSquareCheck)
    _, _, hist = countsAndHistogram(values, check)
    hist
end

function isDistribution(values :: Vector, check :: ChiSquareCheck)
    values_count, intervals_count, hist = countsAndHistogram(values, check)

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
