module Variance
using ValidatedNumerics;
using Distributions;
using Dierckx;

function byConfidenceProbability(average, values, alpha, length)
    gamma_down, gamma_up = getGammas(alpha)
    quantile_down = getChiSquareQuantile(gamma_down, length)
    quantile_up = getChiSquareQuantile(gamma_up, length)

    reduced = mapreduce((x) -> (x - average)^2, +, values)

    @interval(reduced/quantile_down, reduced/quantile_up)
end

function byMeanAbsoluteDeviation(average, values, alpha, length)
    quantile_down, quantile_up = getMeanAbsDeviationQuantiles(alpha, length)

    mean_abs_deviation = mapreduce((x) -> abs(x - average), +, values)/length

    @interval(mean_abs_deviation/quantile_down, mean_abs_deviation/quantile_up)
end

function byPointVariance(variance, alpha, length)
    gamma_down, gamma_up = getGammas(alpha)
    quantile_down = getChiSquareQuantile(gamma_down, length)
    quantile_up = getChiSquareQuantile(gamma_up, length)

    fixed_variance = sqrt(variance) * (1 + 0.254 / (length - 1))
    term = (length - 1) * fixed_variance*fixed_variance

    @interval(term / quantile_down, term / quantile_up)
end

function getGammas(alpha)
    (1 + alpha) * 0.5, (1 - alpha) * 0.5
end

function getChiSquareQuantile(value, length)
    d = Chisq(length - 1)
    quantile(d, value)
end

function getMeanAbsDeviationQuantiles(alpha, length)
    if alpha == 0.90
        knots = [
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
    elseif alpha == 0.95
        knots = [
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
    elseif alpha == 0.99
        knots = [
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
        error("Value=$alpha is not allowed")
    end
    x_knots = knots[:, 1]
    if length > last(x_knots)
        error("Error $length too high")
    end

    y_down_knots = knots[:, 2]
    down_itp = Spline1D(x_knots, y_down_knots, k=1, bc="extrapolate")

    y_up_knots = knots[:, 3]
    up_itp = Spline1D(x_knots, y_up_knots, k=1, bc="extrapolate")
    evaluate(down_itp, length), evaluate(up_itp, length)
end

end
