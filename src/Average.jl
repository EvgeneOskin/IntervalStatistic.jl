module Average
using ValidatedNumerics;
using Distributions;
using Dierckx;

function byKnownVariance(average, variance, alpha, length)
    gamma = getGamma(alpha)
    quantile_gamma = normalQuantile(gamma)
    delta = quantile_gamma * sqrt(variance) / sqrt(length)
    average + @interval(-delta, delta)
end

function byUnknownVariance(average, values, alpha, length)
    gamma = getGamma(alpha)
    variance = mapreduce((x) -> (x - average)^2, +, values) / (length - 1)
    quantile_gamma = studentQuatile(gamma, length - 1)
    delta = quantile_gamma * sqrt(variance) / sqrt(length)
    average + @interval(-delta, delta)
end

function byMeanAbsDeviation(average, values, alpha, length)
    gamma = getGamma(alpha)
    mean_abs_deviation = mapreduce((x) -> abs(x - average), +, values) / length
    quantile_gamma = getMeanAbsDeviationQuantile(gamma, length)
    delta = mean_abs_deviation * quantile_gamma
    average + @interval(-delta, delta)
end

function byInterQuartileWidth(values, alpha, length)
    gamma = getGamma(alpha)
    ordered_values = sort(values)
    interquartile_width = getInterQuartileWidth(ordered_values, length)
    interquertile_quitile = getInterQuartileQuintile(gamma, length)
    mediam = getMediamOfSorted(ordered_values, length)
    delta = interquartile_width * interquertile_quitile
    mediam + @interval(-delta, delta)
end

function getGamma(alpha)
    (1.0 + alpha) * 0.5
end

function getInterQuartileWidth(values, length)
    quartile_index = div((length + 1), 4)
    values[3*quartile_index] - values[quartile_index]
end

function getMediamOfSorted(values, length)
    if isodd(length)
        0.5*(values[div(length, 2)] + values[div(length, 2) + 1])
    else
        values[div((length + 1), 2)]
    end
end

function getInterQuartileQuintile(value, length)
    knots = [
    11 0.470 0.623 0.876
    15 0.400 0.514 0.678
    19 0.354 0.448 0.573
    23 0.321 0.402 0.506
    27 0.296 0.367 0.458
    31 0.276 0.341 0.422
    35 0.260 0.319 0.393
    39 0.246 0.301 0.369
    43 0.234 0.289 0.349
    47 0.224 0.273 0.332
    51 0.215 0.261 0.317
    ]
    x_knots = knots[:, 1]
    if value == 0.95
        y_knots = knots[:, 2]
    elseif value == 0.975
        y_knots = knots[:, 3]
    elseif value == 0.99
        y_knots = knots[:, 4]
    else
        error("Value=$value is not allowed")
    end
    itp = Spline1D(x_knots, y_knots, k=1, bc="extrapolate")
    evaluate(itp, length)
end

function getMeanAbsDeviationQuantile(value, length)
    if value != 0.95 && value != 0.975
        error("Value=$value is not allowed")
    end
    knots = [
    2 12.71;
    3 3.45;
    4 2.16;
    5 1.66;
    6 1.40;
    7 1.21;
    8 1.09;
    9 1.00;
    10 0.93;
    11 0.87;
    12 0.82;
    13 0.78;
    14 0.75;
    20 0.71;
    25 0.60;
    30 0.48;
    40 0.41;
    60 0.33;
    120 0.23;
    ]
    x_knots = knots[:, 1]
    y_knots = knots[:, 2]
    itp = Spline1D(x_knots, y_knots, k=1, bc="extrapolate")
    evaluate(itp, length)
end

function studentQuatile(value, freedomDegreeNumber)
    d = TDist(freedomDegreeNumber)
    quantile(d, value)
end

function normalQuantile(value)
    d = Normal()
    quantile(d, value)
end
end
