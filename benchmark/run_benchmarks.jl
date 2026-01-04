using Pkg
Pkg.activate(temp=true)
Pkg.develop(path=".")
Pkg.add("BenchmarkTools")

using Hurst
using BenchmarkTools
using Random

Random.seed!(1234)

println("Setting up benchmark...")
n = 100_000
X = accumulate(+, randn(n))
τ_range = 1:100
q_range = 0.0:0.1:2.0

println("Benchmarking generalised_hurst_range with use_threading=true:")
b_threaded = @benchmark generalised_hurst_range($X, $τ_range, $q_range, use_threading=true)
display(b_threaded)

println("\nBenchmarking generalised_hurst_range with use_threading=false:")
b_serial = @benchmark generalised_hurst_range($X, $τ_range, $q_range, use_threading=false)
display(b_serial)
