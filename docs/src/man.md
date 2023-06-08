# Using Hurst.jl

`Hurst.jl` is designed to be used in two ways: directly calculating exponents for individual time series, or to systematically analyse large number of time series.

## Individual series 
To calculate an individual series, just call `hurst_exponent` on your time series data along with the range of ``\tau`` that you would like to calculate over: 

```@example
using Hurst
X = accumulate(+, randn(1000))
tau_range = 1:19
hurst_exponent(X, tau_range)
```

This always returns the standard error of the Hurst statistic too.

Given that Hurst exponents are one measure of a series' scaling properties, Generalised Hurst Exponents (GHEs) can be used to detect more complex scaling properties and infer multiscaling behaviour. For example, to calculate the Hurst exponent over a range of moments `q_range`:

```@example
using Hurst
q_range = 0.1:0.1:1
tau_range = 30:250
X = accumulate(+, randn(1000))
generalised_hurst_range(X, tau_range, q_range)
```

## Many series 
For larger investigations, we can usefully reuse most of the memory used in the hurst analysis. For example, suppose we want to calculate the ``\zeta(q)`` estimator for `A` and `B`. To calculate the exponent for `A`, `Hurst.jl` will create dependent and independet variables for a regression calculation, that are normally discarded after use...

