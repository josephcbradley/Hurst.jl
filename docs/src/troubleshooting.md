# Troubleshooting 
The most likely problem to encounter when using this package is ```NaN``` outputs from ```hurst_exponent(...)```. This normally happens when the data you have provided is not self-similar in the way the package expects. Before filing an issue on Github, please read the advice below. 

## TL;DR:
Do the following with your data ```X```:

1. Calculate ```τ_range = 30:250; Q = [Hurst.qth_abs_moment(X, τ, 1) for τ in τ_range]```
2. Check ```plot(log.(τ_range), log.(Q))```

If the plot isn't roughly linear, you have a problem and you will not get good estimates for $H$.

## Why am I getting weird results for $H$?
The easiest way to demonstrate what is probably going wrong is with Gaussian data, so we will generate some Gaussian increments ($X(t) \sim \mathcal{N}(0, 1)$) and a Gaussian walk ($Z(t) = \sum_\tau X(\tau), X(\tau) \sim \mathcal{N}(0, 1)$).

```@example troubleshooting_1
using Random 
Random.seed!(123)
using Plots
gaussian_increments = randn(3000)
gaussian_walk = cumsum(gaussian_increments)
plot(
    (
        plot(gaussian_increments, label = "Gaussian increments"),
        plot(gaussian_walk, label = "Gaussian walk")
    )..., layout = (2, 1)
)
```

## Verify that the quantity is scaling 
```Hurst.jl``` assumes your data $X$ is self-similar, by which I mean:

$\left\{X\left(c t_1\right), \ldots, X\left(c t_k\right)\right\} \stackrel{\mathrm{d}}{=}\left\{c^H X\left(t_1\right), \ldots, c^H X\left(t_k\right)\right\}$

for some $H > 0$ and a process $X(t)$. This definition is probably only intuitive to people already familiar with the subject, so one way to visualise this as follows:

1. Define a time range, ```t_range```, a scale ```c``` - let's say ```0.5``` - and a scaling exponent ```H``` (also ```0.5```).
2. Plot the distribution of ```X[c .* t_range]```, i.e. generate lots of ```t_range``` and plot ```X[c .* t_range]``` each time.
3. Similarly plot the distribution of ```c^H .* X[t_range]```
4. Compare the distributions.

```@example troubleshooting_1 
using Hurst

c = 0.1
H = 0.5
ntrials = 200
alpha = 0.2
starting_t = rand(100:200)
t_range = collect(range(starting_t, length = 20, step = 50))
scaled_t_range = round.(Int64, c .* t_range)
lhs_plt = plot(gaussian_walk[scaled_t_range], label = "ctₖ", color = :blue, alpha = alpha)
rhs_plt = plot(c^H .* gaussian_walk[t_range], label = "c^H", color = :blue, alpha = alpha)
for n = 1:ntrials - 1
    starting_t = rand(100:500)
    t_range .= collect(range(starting_t, length = 20, step = 100))
    scaled_t_range .= round.(Int64, c .* t_range)
    plot!(lhs_plt, gaussian_walk[scaled_t_range], label = nothing, color = :blue, alpha = alpha)
    plot!(rhs_plt, c^H .* gaussian_walk[t_range], label = nothing, color = :blue, alpha = alpha)
end		
plot((lhs_plt, rhs_plt)..., layout = (2, 1))
```

These plots are not perfectly similar - we're only using 200 samples in this case - but hopefully this scaling property is clear enough. When we scale time by $c$, if we scale our output by $c^H$ we obtain a process that is statistically similar.

Now let's try with something that does not have this nice scaling property, e.g. $\sin(x_i)$ for ```i = 1:3000```:

```@example troubleshooting_1 
X = sin.(1:3000)
starting_t = rand(100:500)
t_range = collect(range(starting_t, length = 20, step = 100))
scaled_t_range = round.(Int64, c .* t_range)
lhs_plt = plot(X[scaled_t_range], label = "ctₖ", color = :blue, alpha = alpha)
rhs_plt = plot(c^H .* X[t_range], label = "c^H", color = :blue, alpha = alpha)
for n = 1:ntrials - 1
    starting_t = rand(100:500)
    t_range .= collect(range(starting_t, length = 20, step = 100))
    scaled_t_range .= round.(Int64, c .* t_range)
    plot!(lhs_plt, X[scaled_t_range], label = nothing, color = :blue, alpha = alpha)
    plot!(rhs_plt, c^H* X[t_range], label = nothing, color = :blue, alpha = alpha)
end		
plot((lhs_plt, rhs_plt)..., layout = (2, 1))
```
Now the distributions look very different! If we try to calculate the Hurst exponent of $\sin(x)$ we will get very strange results. Under the hood, ```Hurst.jl``` considers the moment of the absolute $\tau$-increments of the time series. Call this $Q$:

$Q = E\left[|X(t+\tau)-X(t)|\right]$

where $\tau$ represents some time delay. We assume that it is equal to a power law in $\tau$:

$Q=K\tau^{H}$

We calculate the LHS for various values of $\tau$ and take logs:

$\ln (Q) = H \ln (\tau)+\ln (K)$

Our Hurst exponent $H$ is just the coefficient of this linear regression. We can do this all using the functions in ```Hurst.jl```:

```@example troubleshooting_1
τ_range = 1:19 #not the best choice but good for visualising 
qth_moment_over_tau = [Hurst.qth_abs_moment(gaussian_walk, τ, 1) for τ ∈ τ_range]
scatter(log.(τ_range), log.(qth_moment_over_tau), title = "Gaussian walk",
    xlabel = "ln(τ)", label = "ln qth absolute moment")
using GLM
gaussian_lm = lm(hcat(ones(length(τ_range)), log.(τ_range)), log.(qth_moment_over_tau))
a, b = coef(gaussian_lm)
plot!(lnτ -> a + b*lnτ, log.(τ_range), label = "Linear fit")
```
It is easy to see that a lienar regression on this will be meaningful and will give a good estimate for $H$. 

However, if we do not have self-similar data, this linear relationship between $\ln (\tau)$ and $\ln (Q)$ will not be present, and we will get invalid values for $H$:

```@example troubleshooting_1
qth_moment_over_tau = [Hurst.qth_abs_moment(gaussian_increments, τ, 1) for τ ∈ τ_range]
scatter(log.(τ_range), log.(qth_moment_over_tau), title = "sin(x)",
    xlabel = "ln(τ)", label = "ln qth absolute moment");
increments_lm = lm(hcat(ones(length(τ_range)), log.(τ_range)), log.(qth_moment_over_tau))
lnK, H = coef(increments_lm)
```

```@example troubleshooting_1
plot!(lnτ -> lnK + H*lnτ, log.(τ_range), label = "Linear fit")
```

A negative Hurst exponent is not possible, and only appears when we try to estimate one for data that is not scaling. If you try to calculate the Hurst exponent on data like this, ```Hurst.jl``` will return ```NaN```:

```@example troubleshooting_1
H, SD = hurst_exponent(gaussian_increments, τ_range)
```

At this point, hopefully you understand the motivation for the routine at the start of this page:

1. Calculate ```Q = [Hurst.qth_abs_moment(X, τ, q) for τ in τ_range]```
2. Check ```plot(log.(τ_range), log.(Q))```

If the plot isn't roughly linear, you have a problem and you will not get good estimates for $H$.


```@example troubleshooting_1
q = 1.
Q = [Hurst.qth_abs_moment(gaussian_increments, τ, q) for τ in τ_range]
plt = scatter(log.(τ_range), log.(Q), label = "data")
nasty_lm = lm(hcat(ones(length(τ_range)), log.(τ_range)), log.(Q))
a, b = coef(nasty_lm)
plot!(plt, lnτ -> a + b*lnτ, log.(τ_range), label = "Linear fit")
plt
#not linear!
```
```@example troubleshooting_1
Q = [Hurst.qth_abs_moment(gaussian_walk, τ, q) for τ in τ_range]
plt = scatter(log.(τ_range), log.(Q), label = "data")
nice_lm = lm(hcat(ones(length(τ_range)), log.(τ_range)), log.(Q))
a, b = coef(nice_lm)
plot!(plt, lnτ -> a + b*lnτ, log.(τ_range), label = "Linear fit")
plt
#much better
```