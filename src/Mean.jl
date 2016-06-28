module Mean
using ValidatedNumerics
using Distributions
using Dierckx

abstract BaseEstimator

function estimate(values, method :: BaseEstimator)
    error("Unsupported method")
end

immutable byKnownVariance <: BaseEstimator
    confidence_probability :: Real
    variance :: Real
end

function estimate(values, method::byKnownVariance)
    interval_confidence = getIntervalConfidenceProbability(method.confidence_probability)
    average = mean(values)
    count = length(values)
    quantile = normalQuantile(interval_confidence)
    delta = quantile * sqrt(method.variance) / sqrt(count)
    average + @interval(-delta, delta)
end

immutable byUnknownVariance <: BaseEstimator
    confidence_probability :: Real
end

function estimate(values, method::byUnknownVariance)
    average = mean(values)
    count = length(values)
    interval_confidence = getIntervalConfidenceProbability(method.confidence_probability)
    variance = mapreduce((x) -> (x - average)^2, +, values) / (count - 1)
    quantile = studentQuatile(interval_confidence, count - 1)
    delta = quantile * sqrt(variance) / sqrt(count)
    average + @interval(-delta, delta)
end

immutable byMeanAbsDeviation <: BaseEstimator
    confidence_probability :: Real
end

function estimate(values, method :: byMeanAbsDeviation)
    average = mean(values)
    count = length(values)
    interval_confidence = getIntervalConfidenceProbability(method.confidence_probability)
    mean_abs_deviation = mapreduce((x) -> abs(x - average), +, values) / count
    quantile = getMeanAbsDeviationQuantile(interval_confidence, count)
    delta = mean_abs_deviation * quantile
    average + @interval(-delta, delta)
end

immutable byInterQuartileWidth <: BaseEstimator
    confidence_probability :: Real
end

function estimate(values, method :: byInterQuartileWidth)
    interval_confidence = getIntervalConfidenceProbability(method.confidence_probability)
    count = length(values)
    ordered_values = sort(values)
    interquartile_width = getInterQuartileWidth(ordered_values, count)
    interquertile_quitile = getInterQuartileQuintile(interval_confidence, count)
    mediam = getMediamOfSorted(ordered_values, count)
    delta = interquartile_width * interquertile_quitile
    mediam + @interval(-delta, delta)
end

function getIntervalConfidenceProbability(confidence_probability)
    (1.0 + confidence_probability) * 0.5
end

function getInterQuartileWidth(values, count)
    # TODO this can be replaced with StatsBase.iqr
    quartile_index = div((count + 1), 4)
    values[3*quartile_index] - values[quartile_index]
end

function getMediamOfSorted(values, count)
    if isodd(count)
        0.5*(values[div(count, 2)] + values[div(count, 2) + 1])
    else
        values[div((count + 1), 2)]
    end
end

function getInterQuartileQuintile(value, count)
    const knots = [
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
    evaluate(itp, count)
end

function getMeanAbsDeviationQuantile(value, count)
    if value != 0.95 && value != 0.975
        error("Value=$value is not allowed")
    end
    const knots = [
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
    evaluate(itp, count)
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
