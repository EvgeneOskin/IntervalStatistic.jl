module CheckTests
using Base.Test
using FactCheck
using ValidatedNumerics
using IntervalStatistic
using Distributions

srand(10)

in_3_sigma_interval(mu, sigma) = (i) -> (mu in i) && (i in mu + 3*@interval(-sigma, sigma))

facts("test goodness-fit of standard (0, 1) distribution with true distribution") do
    d = Normal()
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

    context("by chi_square check with Sturges k") do
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.SturgesChiSquareCheck(0.05, d)
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Scott k") do
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.ScottChiSquareCheck(0.05, d)
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Taylor k") do
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.TaylorChiSquareCheck(0.05, d)
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with FreedmanDiaconis k") do
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.FreedmanDiaconisChiSquareCheck(0.05, d)
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Doane k") do
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.DoaneChiSquareCheck(0.05, d)
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Wichard k") do
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.WichardChiSquareCheck(0.05, d)
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square modified for normal dist check with 10 bins") do
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.ChiSquareNormalCheck(0.05, 10, mu, sigma)
        )
        println(result)
        @fact result --> true
    end

    context("by chi_square modified for normal dist check with 11 bins") do
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.ChiSquareNormalCheck(0.05, 11, mu, sigma)
        )
        println(result)
        @pending @fact result --> true
    end
end

facts("test goodness-fit of standard (0, 1) distribution") do
    d = Normal()
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

    context("by chi_square check with Sturges k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.SturgesChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Scott k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.ScottChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Taylor k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.TaylorChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with FreedmanDiaconis k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.FreedmanDiaconisChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Doane k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.DoaneChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Wichard k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.WichardChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square modified for normal dist check with 10 bins") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.ChiSquareNormalCheck(0.05, 10, mu, sigma)
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square modified for normal dist check with 11 bins") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.ChiSquareNormalCheck(0.05, 11, mu, sigma)
        )
        println(result)
        @pending @fact result --> true
    end
end

facts("test goodness-fit of standard (0, 1) distribution with samlping=500") do
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
        @pending @fact result --> true
    end
end

facts("test goodness-fit of normal mu=3, sigma=0.1 distribution") do
    d = Normal(3, 0.1)
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

    context("by chi_square check with Dahiya k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.SturgesChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Scott k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.ScottChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Taylor k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.TaylorChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with FreedmanDiaconis k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.FreedmanDiaconisChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Doane k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.DoaneChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Wichard k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.WichardChiSquareCheck(0.05, Normal(mu, sigma))
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square modified for normal dist check with 10 bins") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.ChiSquareNormalCheck(0.05, 10, mu, sigma)
        )
        println(result)
        @fact result --> true
    end

    context("by chi_square modified for normal dist check with 11 bins") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.ChiSquareNormalCheck(0.05, 11, mu, sigma)
        )
        println(result)
        @fact result --> true
    end
end


facts("test goodness-fit of normal mu=3, sigma=0.1 distribution with true distribution") do
    d = Normal(3, 0.1)
    length = 100
    values = rand(d, length)
    mu, sigma = params(d)
    average = mean(values)

    context("by chi_square check with Dahiya k") do
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.SturgesChiSquareCheck(0.05, d)
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Scott k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.ScottChiSquareCheck(0.05, d)
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Taylor k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.TaylorChiSquareCheck(0.05, d)
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with FreedmanDiaconis k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.FreedmanDiaconisChiSquareCheck(0.05, d)
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Doane k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.DoaneChiSquareCheck(0.05, d)
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square check with Wichard k") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.WichardChiSquareCheck(0.05, d)
        )
        println(result)
        @pending @fact result --> true
    end

    context("by chi_square modified for normal dist check with 10 bins") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.ChiSquareNormalCheck(0.05, 10, mu, sigma)
        )
        println(result)
        @fact result --> true
    end

    context("by chi_square modified for normal dist check with 11 bins") do
        mu = mean(values)
        sigma = sqrt(var(values))
        result = IntervalStatistic.isDistribution(
            values,
            IntervalStatistic.Check.ChiSquareNormalCheck(0.05, 11, mu, sigma)
        )
        println(result)
        @fact result --> true
    end
end

facts("test goodness-fit of normal mu=3, sigma=0.1 distribution with sampling=500") do
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
        @pending @fact result --> true
    end
end


end
