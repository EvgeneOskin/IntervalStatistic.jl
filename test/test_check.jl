module CheckTests
using Base.Test
using FactCheck
using ValidatedNumerics
using IntervalStatistic
using Distributions

srand(10)

in_3_sigma_interval(mu, sigma) = (i) -> (mu in i) && (i in mu + 3*@interval(-sigma, sigma))

facts("estimate average of standard (0, 1) distribution") do
    d = Normal()
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

    context("by simple chi_square check") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.SimpleChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @pending @fact result --> true
    end
    context("by chi_square check with Dahiya k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.DahiyaChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @fact result --> true
    end
end

facts("estimate average of standard (0, 1) distribution with samlping=500") do
    d = Normal()
    length = 500
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

    context("by chi_square check with k for Large n") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.LargeNChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @fact result --> true
    end
end

facts("estimate average of normal mu=3, sigma=0.1 distribution") do
    d = Normal(3, 0.1)
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

    context("by simple chi_square check") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values, IntervalStatistic.Check.SimpleChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @fact result --> true
    end
    context("by chi_square check with Dahiya k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.DahiyaChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @fact result --> true
    end
end

facts("estimate average of normal mu=3, sigma=0.1 distribution with sampling=500") do
    d = Normal(3, 0.1)
    length = 500
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

    context("by chi_square check with k for Large n") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.LargeNChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @fact result --> true
    end
end


end
