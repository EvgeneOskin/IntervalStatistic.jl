module CheckTests
using Base.Test;
using FactCheck;
using ValidatedNumerics;
using IntervalStatistic;
using Distributions;

srand(10)

in_3_sigma_interval(mu, sigma) = (i) -> (mu in i) && (i in mu + 3*@interval(-sigma, sigma))

facts("estimate average of standard (0, 1) distribution") do
    d = Normal()
    length = 1000
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

    context("by simple chi_square check") do
        result = IntervalStatistic.isDistribution(
            values, IntervalStatistic.Check.SimpleChiSquareCheck(0.10)
        )
        println(result)
        @fact result --> true
    end
end

facts("estimate average of standard (0, 1) distribution") do
    d = Normal(3, 0.1)
    length = 1000
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

    context("by simple chi_square check") do
        result = IntervalStatistic.isDistribution(
            values, IntervalStatistic.Check.SimpleChiSquareCheck(0.10)
        )
        println(result)
        @fact result --> true
    end
end

end

