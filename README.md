# Interval Statistic

[![Build Status](https://travis-ci.org/EvgeneOskin/IntervalStatistic.jl.svg?branch=feature%2Ftune-up-travis)](https://travis-ci.org/EvgeneOskin/IntervalStatistic.jl)

# Intallation

```julia
Pkg.add("IntervalStatistic")
```

# Usage
```julia
julia> using IntervalStatistic
julia> using Distributions
julia> d = Normal()
julia> length = 10
julia> values = rand(d, length)
julia> mu, sigma = params(d)
julia> average = reduce(+, values) / length
julia> confidence_probability = 0.95
julia> # Average estimations
julia> IntervalStatistic.Average.byKnownVariance(average, sigma*sigma, confidence_probability, length)
julia> IntervalStatistic.Average.byInterQuartileWidth(values, confidence_probability, length)
julia> IntervalStatistic.Average.byUnknownVariance(average, values, confidence_probability, length)
julia> IntervalStatistic.Average.byMeanAbsDeviation(average, values, confidence_probability, length)
julia> # Variance estimations
julia> IntervalStatistic.Variance.byConfidenceProbability(average, values, confidence_probability, length)
julia> variance = mapreduce((x) -> (x - average)^2, +, values) / length
julia> IntervalStatistic.Variance.byPointVariance(variance, confidence_probability, length)
julia> IntervalStatistic.Variance.byMeanAbsoluteDeviation(average, values, confidence_probability, length)

```

Julia package to keep functions to estimate interval average and variance.
