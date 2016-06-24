module AverageTests
using Base.Test;
using FactCheck;
using ValidatedNumerics;
using IntervalStatistic;
using Distributions;

facts("By known Variance standard (0, 1) distribution") do
    d = Normal()
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length
    result = IntervalStatistic.Average.byKnownVariance(average, sigma*sigma, 0.95, length)
    println(result)
    @fact mu --> less_than(result.hi)
    @fact mu --> greater_than(result.lo)
    @fact mu + 3 * sigma --> greater_than(result.hi)
    @fact mu - 3 * sigma --> less_than(result.lo)
end

facts("By known Variance mu=3 sigma=0.1 distribution") do
    d = Normal(3, 0.1)
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length
    result = IntervalStatistic.Average.byKnownVariance(average, sigma*sigma, 0.95, length)
    println(result)
    @fact mu --> less_than(result.hi)
    @fact mu --> greater_than(result.lo)
    @fact mu + 3 * sigma --> greater_than(result.hi)
    @fact mu - 3 * sigma --> less_than(result.lo)
end

facts("By unknown Variance standard (0, 1) distribution") do
    d = Normal()
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length
    result = IntervalStatistic.Average.byUnknownVariance(average, values, 0.95, length)
    println(result)
    @fact mu --> less_than(result.hi)
    @fact mu --> greater_than(result.lo)
    @fact mu + 3 * sigma --> greater_than(result.hi)
    @fact mu - 3 * sigma --> less_than(result.lo)
end

facts("By unknown Variance mu=3 sigma=0.1 distribution") do
    d = Normal(3, 0.1)
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = reduce(+, values) / length
    result = IntervalStatistic.Average.byUnknownVariance(average, values, 0.95, length)
    println(result)
    @fact mu --> less_than(result.hi)
    @fact mu --> greater_than(result.lo)
    @fact mu + 3 * sigma --> greater_than(result.hi)
    @fact mu - 3 * sigma --> less_than(result.lo)
end

end

