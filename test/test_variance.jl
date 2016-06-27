module VarianceTests
using Base.Test;
using FactCheck;
using ValidatedNumerics;
using IntervalStatistic;
using Distributions;

srand(10)

contain(value) = (i) -> value in i

facts("estimate variance of standard (0, 1) distribution") do
    d = Normal()
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

    context("by confidence probability") do
        result = IntervalStatistic.var(values, IntervalStatistic.Variance.byConfidenceProbability(0.95))
        println(result)
        @fact result --> contain(sigma^2)
    end

    context("by point variance") do
        variance = var(values)
        result = IntervalStatistic.var(values, IntervalStatistic.Variance.byPointVariance(0.95, variance))
        println(result)
        @fact result --> contain(variance)
    end
end

facts("estimate variance of normanl mu=3 sigma=0.1 distribution") do
    d = Normal(3, 0.1)
    length = 10
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

    context("by confidence probability") do
        result = IntervalStatistic.var(values, IntervalStatistic.Variance.byConfidenceProbability(0.95))
        println(result)
        @fact result --> contain(sigma^2)
    end

    context("by point variance") do
        variance = var(values)
        result = IntervalStatistic.var(values, IntervalStatistic.Variance.byPointVariance(0.95, variance))
        println(result)
        @fact result --> contain(variance)
    end
end

facts("estimate variance of standard (0, 1) distribution with 10 samples") do
    d = Normal()
    length = 10
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

    context("by mean abs deviation") do
        result = IntervalStatistic.var(values, IntervalStatistic.Variance.byMeanAbsoluteDeviation(0.95))
        println(result)
        @fact result --> contain(sigma^2)
    end
end

facts("estimate variance of normanl mu=3 sigma=0.1 distribution with 10 samples") do
    d = Normal(3, 0.1)
    length = 10
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

   context("by mean abs deviation") do
        result = IntervalStatistic.var(values, IntervalStatistic.Variance.byMeanAbsoluteDeviation(0.95))
        println(result)
        @fact result --> contain(sigma^2)
    end
end

end
