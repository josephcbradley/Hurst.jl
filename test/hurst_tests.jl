@testset "Hurst exponent tests" begin 
    gaussian_walk = accumulate(+, randn(10000))
    H_g, H_s = hurst_exponent(gaussian_walk, 1:50)
    @test isapprox(H_g, 0.5, atol = 0.05) #this test will always be rough!

    # test method error
    g_matrix = reshape(gaussian_walk, (100, 100))
    @test_throws MethodError hurst_exponent(g_matrix, 1:50)
end

@testset "qth-moment tests" begin 
    A = collect(1:1:100)
    @test Hurst.qth_abs_moment(A, 1, 1) == 1.
    @test Hurst.qth_abs_moment(A, 2, 1) == 2.
    @test_throws DimensionMismatch Hurst.qth_abs_moment(A, 200, 2)
end