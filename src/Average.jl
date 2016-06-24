module Average
using ValidatedNumerics;
using Distributions;
using Dierckx;

function byKnownVariance(average, variance, alpha, length)
    gamma = getGamma(alpha)
    quantile_gamma = normalQuantile(gamma)
    delta = quantile_gamma * sqrt(variance) / sqrt(length)
    @interval(average - delta, average + delta)
end

function byUnknownVariance(average, values, alpha, length)
    gamma = getGamma(alpha)
    variance = mapreduce((x) -> (x - average)^2, +, values) / (length - 1)
    quantile_gamma = studentQuatile(gamma, length - 1)
    delta = quantile_gamma * sqrt(variance) / sqrt(length)
    @interval(average - delta, average + delta)
end

function byMeanAbsDeviation(average, values, alpha, length)
    gamma = getGamma(alpha)
    mean_abs_deviation = mapreduce((x) -> abs(x - average), +, values) / length
    quantile_gamma = getMeanAbsDeviationQuantile(gamma, length)
    delta = average_abs_deviation * quantile_gamma
    @interval(average - delta, average + delta)
end

function getGamma(alpha)
    (1.0 + alpha) * 0.5
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
