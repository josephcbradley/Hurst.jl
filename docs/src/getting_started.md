# Getting Started
`Hurst.jl` allows you to calculate Generalised Hurst Exponents (GHEs) in a robust and performant way.

## Install Hurst.jl
Before anything else, install `Hurst.jl` in the usual way:

```@example getting_started_1 
using Pkg
Pkg.add("Hurst")
```

## Get some data
The easiest way to demonstrate how to use ```Hurst.jl``` is with Gaussian data, so we will generate a Gaussian walk ($Z(t) = \sum_\tau X(\tau), X(\tau) \sim \mathcal{N}(0, 1)$).

```@example getting_started_1
using Hurst, Plots, Random
Random.seed!(1)  
X = cumsum(randn(3000))
plot(X, label = "Gaussian walk", xlabel = "t", ylabel = "X(t)")
```

## Define your $\tau$ range 
```Hurst.jl``` calculates the Hurst exponent $H$ by estimating the moments of the absolute increments over some time delay $\tau$. If that sounds like gibberish, take a look at the [Troubleshooting](@ref) section which explains things in more detail. ```τ_range = 30:250``` is a good place to start:

```@example getting_started_1
τ_range = 30:250
```
## Calculate $H$
Calclating $H$ is as simple as:

```@example getting_started_1
H, SD = hurst_exponent(X, τ_range)
```

Notice that ```hurst_exponent``` returns a ```1x2``` array. This is because we want to be able to calculate the Generalised Hurst Exponent (GHE) quickly as well and have a consistent API:


```@example getting_started_1
q_range = 0.1:0.1:1.
GHE_data = generalised_hurst_range(X, τ_range, q_range)
```

```@example getting_started_1
plot(q_range, GHE_data[:,1], yerror = GHE_data[:, 2], xlabel = "q",
    ylabel = "H(q)", label = nothing, ylims = (0, 1))
```

## Further analysis 
If you are trying to reproduce some results from papers that use the GHE, you might be interested in the quantity $\zeta(q) = qH(q)$:

```@example getting_started_1
zeta_data = Hurst.zeta_estimator_range(X, τ_range, q_range)
plot(q_range, zeta_data[:,1], yerror = zeta_data[:, 2], xlabel = "q",
    ylabel = "ζ(q)", label = nothing)
```

If you want to see if the data is multi-scaling, you can fit parabola to this data with ```Polynomials.jl```:

```@example getting_started_1
using Polynomials
quadratic_fit = fit(q_range, zeta_data[:, 1], 2)
scatter!(q_range, zeta_data[:, 1], label = nothing)
plot!(quadratic_fit, extrema(q_range)..., label = "Quadratic fit")
```

```@example getting_started_1
#first coefficient should be almost zero...
coeffs(quadratic_fit)
```
