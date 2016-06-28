module MeanTests
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
    average = mean(values)

    context("by known variance") do
        result = IntervalStatistic.mean(values, IntervalStatistic.Mean.byKnownVariance(0.95, sigma*sigma))
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end

    context("by interquartile width") do
        result = IntervalStatistic.mean(values, IntervalStatistic.Mean.byInterQuartileWidth(0.95))
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end

    context("by unknown variance") do
        result = IntervalStatistic.mean(values, IntervalStatistic.Mean.byUnknownVariance(0.95))
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end

    context("by mean abs deviation") do
        result = IntervalStatistic.mean(values, IntervalStatistic.Mean.byMeanAbsDeviation(0.95))
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end
end

facts("estimate average of normal mu=3, sigma=0.1 distribution") do
    d = Normal(3, 0.1)
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

    context("by known variance") do
        result = IntervalStatistic.mean(values, IntervalStatistic.Mean.byKnownVariance(0.95, sigma*sigma))
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end

    context("by interquartile width") do
        result = IntervalStatistic.mean(values, IntervalStatistic.Mean.byInterQuartileWidth(0.95))
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end

    context("by unknown variance") do
        result = IntervalStatistic.mean(values, IntervalStatistic.Mean.byUnknownVariance(0.95))
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end

    context("by mean abs deviation") do
        result = IntervalStatistic.mean(values, IntervalStatistic.Mean.byMeanAbsDeviation(0.95))
        println(result)
        @fact result --> in_3_sigma_interval(mu, sigma)
    end
end

end
