module Variance
using ValidatedNumerics;
using Distributions;
using Dierckx;

abstract BaseEstimator

function estimate(values, method::BaseEstimator)
    error("Unsupported method $method")
end

immutable byConfidenceProbability <: BaseEstimator
    confidence_probability :: Real
end

function estimate(values, method :: byConfidenceProbability)
    average = mean(values)
    count = length(values)
    gamma_down, gamma_up = getGammas(method.confidence_probability)
    quantile_down = getChiSquareQuantile(gamma_down, count)
    quantile_up = getChiSquareQuantile(gamma_up, count)

    reduced = mapreduce((x) -> (x - average)^2, +, values)

    @interval(reduced/quantile_down, reduced/quantile_up)
end

immutable byMeanAbsoluteDeviation <: BaseEstimator
    confidence_probability :: Real
end

function estimate(values, method :: byMeanAbsoluteDeviation)
    average = mean(values)
    count = length(values)
    quantile_down, quantile_up = getMeanAbsDeviationQuantiles(method.confidence_probability, count)

    mean_abs_deviation = mapreduce((x) -> abs(x - average), +, values)/count

    @interval(mean_abs_deviation/quantile_down, mean_abs_deviation/quantile_up)^2
end

immutable byPointVariance <: BaseEstimator
    confidence_probability :: Real
    variance :: Real
end

function estimate(values, method::byPointVariance)
    count = length(values)
    gamma_down, gamma_up = getGammas(method.confidence_probability)
    quantile_down = getChiSquareQuantile(gamma_down, count)
    quantile_up = getChiSquareQuantile(gamma_up, count)

    fixed_variance = method.variance * (1 + 0.254 / (count - 1))^2
    term = (count - 1) * fixed_variance

    @interval(term / quantile_down, term / quantile_up)
end

function getGammas(confidence_probability)
    (1 + confidence_probability) * 0.5, (1 - confidence_probability) * 0.5
end

function getChiSquareQuantile(value, count)
    d = Chisq(count - 1)
    quantile(d, value)
end

function getMeanAbsDeviationQuantiles(confidence_probability, count)
    if confidence_probability == 0.90
        const knots = [
        2 1.386 0.044;
        3 1.276 0.166;
        4 1.224 0.254;
        5 1.187 0.305;
        6 1.158 0.360;
        7 1.135 0.394;
        8 1.116 0.422;
        9 1.100 0.445;
        10 1.086 0.464;
        ]
    elseif confidence_probability == 0.95
        const knots = [
        2 1.585 0.022;
        3 1.417 0.116;
        4 1.344 0.199;
        5 1.292 0.260;
        6 1.253 0.306;
        7 1.222 0.342;
        8 1.196 0.372;
        9 1.175 0.396;
        10 1.156 0.417;
        ]
    elseif confidence_probability == 0.99
        const knots = [
        2 1.985 0.004;
        3 1.703 0.073;
        4 1.590 0.145;
        5 1.507 0.203;
        6 1.445 0.250;
        7 1.397 0.287;
        8 1.358 0.318;
        9 1.326 0.344;
        10 1.299 0.366;
        ]
    else
        error("Value=$confidence_probability is not allowed")
    end
    x_knots = knots[:, 1]
    if count > last(x_knots)
        error("Error $count too high")
    end

    y_down_knots = knots[:, 2]
    down_itp = Spline1D(x_knots, y_down_knots, k=1, bc="extrapolate")

    y_up_knots = knots[:, 3]
    up_itp = Spline1D(x_knots, y_up_knots, k=1, bc="extrapolate")
    evaluate(down_itp, count), evaluate(up_itp, count)
end

end
