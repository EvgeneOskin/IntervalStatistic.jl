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

immutable ChiSquareNormalCheck <: ChiSquareCheck
    # Dahiya R.C., Gurland J., How many classes in the Pearson chi-square test? – 1973
    error_probability :: Real
    bin_count :: Int
    distribution
end

function ChiSquareNormalCheck(error_probability :: Real, bin_count :: Int, mu :: Real, sigma :: Real)
    ChiSquareNormalCheck(error_probability, bin_count, Normal(mu, sigma))
end

function fitHistogram(values :: Vector, check :: ChiSquareNormalCheck)
    mu, sigma = params(check.distribution)
    edges = histEdgesForNormalDistribution(check.bin_count, mu, sigma)
    hist = fit(Histogram, values, edges)
end

function histEdgesForNormalDistribution(bin_count, mu, sigma)
    const general_edges_per_bin_count = Dict(
        3=> [-0.4307],
        4=> [-0.6745, 0],
        5=> [-0.8416, -0.2533],
        6=> [-0.9674, -0.4307, 0],
        7=> [-1.0676, -0.5659, -0.1800],
        8=> [-1.1503, -0.6745, -0.3186, 0],
        9=> [-1.2206, -0.7647, -0.4307, -0.1397],
        10=> [-1.2816, -0.8416, -0.5244, -0.2533, 0],
        11=> [-1.3352, -0.9085, -0.6046, -0.3488, -0.1142],
        12=> [-1.3830, -0.9674, -0.6745, -0.4307, -0.2194, 0],
        13=> [-1.4261, -1.0201, -0.7363, -0.5024, -0.2934, -0.0966],
        14=> [-1.4652, -1.0676, -0.7916, -0.5660, -0.3661, -0.1800, 0],
        15=> [-1.5011, -1.1108, -0.8416, -0.6229, -0.4307, -0.2533, -0.0837]
    )

    semi_edges = general_edges_per_bin_count[bin_count]
    semi_edges_length = length(semi_edges)

    edge_symenty = semi_edges_length + (iseven(bin_count) ? 0 : 1)
    calculate_edge = (semi_edges, i) -> begin
        i <= semi_edges_length ? semi_edges[i] : -semi_edges[end - (i - edge_symenty)]
    end

    edges = [calculate_edge(semi_edges, i) for i in 1:(bin_count - 1)]
    far_left_edge, far_right_edge = -5, 5
    edges = vcat([far_left_edge], edges, [far_right_edge])
    mu + sigma*edges
end

function chiSquareCriticalValue(bin_count, check :: ChiSquareNormalCheck)
    critical_values_per_bin_count_and_probability = Dict(
        (3, 0.10)=> 2.371, (3, 0.05)=> 3.248, (3, 0.01)=> 5.418,
        (4, 0.10)=> 3.928, (4, 0.05)=> 5.107, (4, 0.01)=> 7.917,
        (5, 0.10)=> 5.442, (5, 0.05)=> 6.844, (5, 0.01)=> 10.075,
        (6, 0.10)=> 6.905, (6, 0.05)=> 8.479, (6, 0.01)=> 12.021,
        (7, 0.10)=> 8.322, (7, 0.05)=> 10.083, (7, 0.01)=> 13.837,
        (8, 0.10)=> 9.703, (8, 0.05)=> 11.543, (8, 0.01)=> 15.567,
        (9, 0.10)=> 11.055, (9, 0.05)=> 13.007, (9, 0.01)=> 17.234,
        (10, 0.10)=> 12.384, (10, 0.05)=> 14.438, (10, 0.01)=> 18.852,
        (11, 0.10)=> 13.694, (11, 0.05)=> 15.843, (11, 0.01)=> 20.431,
        (12, 0.10)=> 14.688, (12, 0.05)=> 17.226, (12, 0.01)=> 21.977,
        (13, 0.10)=> 16.267, (13, 0.05)=> 19.589, (13, 0.01)=> 23.495,
        (14, 0.10)=> 17.535, (14, 0.05)=> 19.937, (14, 0.01)=> 24.990,
        (15, 0.10)=> 18.792, (15, 0.05)=> 21.270, (15, 0.01)=> 26.464
    )
    bin_coint_and_probability = (bin_count, check.error_probability)
    critical_values_per_bin_count_and_probability[bin_coint_and_probability]
end


function fitHistogram(values :: Vector, check :: ChiSquareCheck)
    values_count = length(values)
    bin_count = calculateBinCount(values, check)
    hist = fit(Histogram, values, nbins=bin_count)
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
    hist = fitHistogram(values, check)
    bin_count = length(hist.weights)
    intervals_and_weights = zip(edgesToIntervals(hist.edges[1]), hist.weights)
    (values_count, bin_count, intervals_and_weights)
end

function fitHistogram(values :: Vector, check :: ChiSquareCheck)
    values_count = length(values)
    bin_count = calculateBinCount(values, check)
    hist = fit(Histogram, values, nbins=bin_count)
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
    values_count, bin_count, hist = countsAndHistogram(values, check)

    chi_square = mapreduce(
        (x) -> chiSquareEstimate(x, values_count, check), +,
        hist
    )

    chi_square_quantile = chiSquareCriticalValue(bin_count, check)
    chi_square <= chi_square_quantile
end

function chiSquareCriticalValue(bin_count, check :: ChiSquareCheck)
    freedom_degree_number = bin_count - 1
    getChiSquareQuantile(check.error_probability, freedom_degree_number)
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
