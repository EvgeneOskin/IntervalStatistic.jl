module VarianceTests
using Base.Test;
using FactCheck;
using ValidatedNumerics;
using IntervalStatistic;
using Distributions;

facts("estimate variance By confidence probability standard (0, 1) distribution") do
    d = Normal()
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length
    result = IntervalStatistic.Variance.byConfidenceProbability(average, values, 0.95, length)
    println(result)
    @fact sigma*sigma --> less_than(result.hi)
    @fact sigma*sigma --> greater_than(result.lo)
end

facts("estimate variance By confidence probability mu=3 sigma=0.1 distribution") do
    d = Normal(3, 0.1)
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length
    result = IntervalStatistic.Variance.byConfidenceProbability(average, values, 0.95, length)
    println(result)
    @fact sigma*sigma --> less_than(result.hi)
    @fact sigma*sigma --> greater_than(result.lo)
end

facts("estimate variance By point variance standard (0, 1) distribution") do
    d = Normal()
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length
    variance = mapreduce((x) -> (x - average)^2, +, values) / length
    result = IntervalStatistic.Variance.byPointVariance(variance, 0.95, length)
    println(result)
    @fact sigma*sigma --> less_than(result.hi)
    @fact sigma*sigma --> greater_than(result.lo)
end

facts("estimate variance By point variance mu=3 sigma=0.1 distribution") do
    d = Normal(3, 0.1)
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length
    variance = mapreduce((x) -> (x - average)^2, +, values) / length
    result = IntervalStatistic.Variance.byPointVariance(variance, 0.95, length)
    println(result)
    @fact sigma*sigma --> less_than(result.hi)
    @fact sigma*sigma --> greater_than(result.lo)
end

facts("estimate variance By mean absolute deviation standard (0, 1) distribution") do
    d = Normal()
    length = 10
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length
    result = IntervalStatistic.Variance.byMeanAbsoluteDeviation(average, values, 0.95, length)
    println(result)
    @fact sigma*sigma --> less_than(result.hi)
    @fact sigma*sigma --> greater_than(result.lo)
end

facts("estimate variance By mean absolute deviation mu=3 sigma=0.1 distribution") do
    d = Normal(3, 0.1)
    length = 10
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length
    result = IntervalStatistic.Variance.byMeanAbsoluteDeviation(average, values, 0.95, length)
    println(result)
    @fact sigma*sigma --> less_than(result.hi)
    @fact sigma*sigma --> greater_than(result.lo)
end

end
