![Hurst.jl logo](https://github.com/josephcbradley/Hurst.jl/blob/main/docs/src/assets/logo.png "Hurst.jl logo")

# Hurst.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://josephcbradley.github.io/Hurst.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://josephcbradley.github.io/Hurst.jl/dev/)
[![Build Status](https://travis-ci.com/josephcbradley/Hurst.jl.svg?branch=main)](https://travis-ci.com/josephcbradley/Hurst.jl)
[![Coverage](https://codecov.io/gh/josephcbradley/Hurst.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/josephcbradley/Hurst.jl)

This package implements methods for estimating [Generalised Hurst Exponents](https://en.wikipedia.org/wiki/Hurst_exponent#Generalized_exponent) (GHEs).

The package is not yet available on the Julia registry, so add it as follows: 

```
using Pkg; Pkg.add("https://github.com/josephcbradley/Hurst.jl")
```

At the moment, the package only implements one method of estimation, and not the most state-of-the-art method at that. In time this will be updated - the goal of the package is to be allow the user to flexibly calculate GHEs in whatever way is most appropriate for the task at hand.

The package is designed to be fast and to facilitate the kind on analysis often seen in the literature. For example, one can calculate the normal Hurst exponent directly:

```
using Hurst
X = accumulate(+, randn(10000));
H = hurst_exponent(X, 1:19)
```

For more methods, including methods to calculate large numbers of Hurst exponents quickly, please see the docs. 

Please note that this package is under development and the interface is likely to change over the coming months. PRs are very welcome.