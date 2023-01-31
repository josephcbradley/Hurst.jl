using Hurst
using Test

@testset "Hurst.jl" begin
    include("linreg_tests.jl")
    include("hurst_tests.jl")
end
