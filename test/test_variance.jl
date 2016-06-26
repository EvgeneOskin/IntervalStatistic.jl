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
    average = reduce(+, values) / length

    context("by confidence probability") do
        result = IntervalStatistic.Variance.byConfidenceProbability(average, values, 0.95, length)
        println(result)
        @fact result --> contain(sigma^2)
    end

    context("by point variance") do
        variance = mapreduce((x) -> (x - average)^2, +, values) / length
        result = IntervalStatistic.Variance.byPointVariance(variance, 0.95, length)
        println(result)
        @fact result --> contain(variance)
    end
end

facts("estimate variance of normanl mu=3 sigma=0.1 distribution") do
    d = Normal(3, 0.1)
    length = 10
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length

    context("by confidence probability") do
        result = IntervalStatistic.Variance.byConfidenceProbability(average, values, 0.95, length)
        println(result)
        @fact result --> contain(sigma^2)
    end

    context("by point variance") do
        variance = mapreduce((x) -> (x - average)^2, +, values) / length
        result = IntervalStatistic.Variance.byPointVariance(variance, 0.95, length)
        println(result)
        @fact result --> contain(variance)
    end
end

facts("estimate variance of standard (0, 1) distribution with 10 samples") do
    d = Normal()
    length = 10
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length

    context("by mean abs deviation") do
        result = IntervalStatistic.Variance.byMeanAbsoluteDeviation(average, values, 0.95, length)
        println(result)
        @fact result --> contain(sigma^2)
    end
end

facts("estimate variance of normanl mu=3 sigma=0.1 distribution with 10 samples") do
    d = Normal(3, 0.1)
    length = 10
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length

   context("by mean abs deviation") do
        result = IntervalStatistic.Variance.byMeanAbsoluteDeviation(average, values, 0.95, length)
        println(result)
        @fact result --> contain(sigma^2)
    end
end

end
