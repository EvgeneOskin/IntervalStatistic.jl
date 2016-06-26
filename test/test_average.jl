module AverageTests
using Base.Test;
using FactCheck;
using ValidatedNumerics;
using IntervalStatistic;
using Distributions;

srand(10)

in_3_sigma_interval(mu, sigma) = (i) -> (mu in i) && (i in mu + 3*@interval(-sigma, sigma))

facts("estimate average of standard (0, 1) distribution") do
    d = Normal()
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length

    context("by known variance") do
        result = IntervalStatistic.Average.byKnownVariance(average, sigma*sigma, 0.95, length)
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end

    context("by interquartile width") do
        result = IntervalStatistic.Average.byInterQuartileWidth(values, 0.95, length)
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end

    context("by unknown variance") do
        result = IntervalStatistic.Average.byUnknownVariance(average, values, 0.95, length)
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end

    context("by mean abs deviation") do
        result = IntervalStatistic.Average.byMeanAbsDeviation(average, values, 0.95, length)
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end
end

facts("estimate average of standard (0, 1) distribution") do
    d = Normal(3, 0.1)
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length

    context("by known variance") do
        result = IntervalStatistic.Average.byKnownVariance(average, sigma*sigma, 0.95, length)
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end

    context("by interquartile width") do
        result = IntervalStatistic.Average.byInterQuartileWidth(values, 0.95, length)
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end

    context("by unknown variance") do
        result = IntervalStatistic.Average.byUnknownVariance(average, values, 0.95, length)
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end

    context("by mean abs deviation") do
        result = IntervalStatistic.Average.byMeanAbsDeviation(average, values, 0.95, length)
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end
end

end
