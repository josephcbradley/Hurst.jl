```@meta
CurrentModule = Hurst
```

# Hurst

This package implements methods for estimating [Generalised Hurst Exponents](https://en.wikipedia.org/wiki/Hurst_exponent#Generalized_exponent) (GHEs).

At the moment, the package only implements one method of estimation, and not the most state-of-the-art method at that. In time this will be updated - the goal of the package is to be allow the user to flexibly calculate GHEs in whatever way is most appropriate for the task at hand.

The package is designed to be fast and to facilitate the kind on analysis often seen in the literature. For example, one can calculate the normal Hurst exponent directly:

```@setup first_example 
using Hurst
X = accumulate(+, randn(10000));
@time H = hurst_exponent(X, 1:19)
```
```@example first_example
using Hurst
X = accumulate(+, randn(10000));
@time H = hurst_exponent(X, 1:19)
```

Or, one can calculate the GHEs for a wide variety of moments:
```@setup second_example 
using Hurst
X = accumulate(+, randn(10000));
```
```@example second_example
using Hurst
tau_range = 1:19
q_range = 0.1:0.1:2.
generalised_hurst_range(X, tau_range, q_range)
```

Hurst exponents (generalised or not) are calculate by performing a regression across a range of values of ``\tau``. It is important to be aware of these values as they can have a significant impact on the results (see [here](https://doi.org/10.1103/PhysRevE.95.042311)). As a result, these values are never supplied by default in functions and must be provided explicitly by the user. Users interest in the details of the calculation are encouraged to look at the above paper and its references.