@testset "Linear regression tests" begin
    #check than linreg works as expected
    f(x) = 2x + 3
    x = collect(-1.0:0.1:1)
    y = f.(x) .+ randn(length(x))
    β_matrix = (x' * x)^-1 * (x' * y)
    β_estimator = Hurst.lsq_estimator(y, x)[1]
    @test β_matrix ≈ β_estimator

    #check that dimension error is thrown 
    @test_throws DimensionMismatch Hurst.lsq_estimator(rand(10), rand(11))
end
