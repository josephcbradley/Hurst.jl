@testset "Hurst exponent tests" begin 
    gaussian_walk = accumulate(+, randn(10000))
    H_g, H_s = hurst_exponent(gaussian_walk, 1:19)
    @test isapprox(H_g, 0.5, atol = 0.1) #this test will always be rough!

    Z_g, Z_s = Hurst.zeta_estimator(gaussian_walk, 1:19, 1)
    @test isapprox(Z_g, 0.5, atol = 0.1) #this test will always be rough!

    z_range = Hurst.zeta_estimator_range(gaussian_walk, 1:19, 0.:0.1:1.)
    @test z_range[1, 1] == 0

    @test generalised_hurst_exponent(gaussian_walk, 1:19, 1)[1] == H_g

    # test method error
    g_matrix = reshape(gaussian_walk, (100, 100))
    @test_throws MethodError hurst_exponent(g_matrix, 1:19)


    #test DimensionMismatch
    Y = rand(10)
    S = rand(11)
    τ_range = 1:12
    @test_throws DimensionMismatch Hurst.zeta_estimator!(Y, S, rand(10), τ_range, 2)


end

@testset "qth-moment tests" begin 
    A = collect(1:1:100)
    @test Hurst.qth_abs_moment(A, 1, 1) == 1.
    @test Hurst.qth_abs_moment(A, 2, 1) == 2.
    @test_throws DimensionMismatch Hurst.qth_abs_moment(A, 200, 2)
end