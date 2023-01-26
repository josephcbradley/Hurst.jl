function zeta_estimator!(Y, S, X, q; τ_range = 1:19)
    #ζ(q) = qH(q) where H is the GHE 
    #let Q = E[|X(t+ τ) - X(t)|^q] = K(q)τ^ζ(q)
    #take logarithms and rearrange
    #ln(Q) = ζ(q)ln(τ) + ln(K(q))
    #to estimate ζ(q), consider the linear regression of
    #Y = β₀ + β₁S where 
    #Y = ln(Q), S = ln(τ), β₀ = ln(K(q)), and β₁ = ζ(q)
    #return β₁ and its std error

    #form Y 
    N = length(τ_range)
    #Y = Vector{Float64}(undef, N)
    #S = Vector{Float64}(undef, N)

    #calculate regression data
    @inbounds for i in 1:N
        τ = τ_range[i]
        Y[i] = log(qth_abs_moment(X, τ, q))
        S[i] = log(τ)
    end

    data = DataFrames.DataFrame(Y = Y, S = S)
    model = GLM.lm(GLM.@formula(Y ~ S), data)
    ζ = GLM.coef(model)[2]
    #get standard error 
    SD = GLM.stderror(model)[2]
    #note ζ(q) = qH(q), so H(q) = ζ(q) / q 
    return ζ, SD
end

function zeta_estimator(X, q, τ_range = 1:19)
    #ζ(q) = qH(q) where H is the GHE 
    #let Q = E[|X(t+ τ) - X(t)|^q] = K(q)τ^ζ(q)
    #take logarithms and rearrange
    #ln(Q) = ζ(q)ln(τ) + ln(K(q))
    #to estimate ζ(q), consider the linear regression of
    #Y = β₀ + β₁S where 
    #Y = ln(Q), S = ln(τ), β₀ = ln(K(q)), and β₁ = ζ(q)
    #return β₁ and its std error

    #form Y 
    N = length(τ_range)
    Y = Vector{Float64}(undef, N)
    S = Vector{Float64}(undef, N)
    zeta_estimator!(Y, S, X, q, τ_range = τ_range)
end

function zeta_estimator_range!(buffer, X; τ_range = 1:19, q_range = 0.:0.1:1.)
    L = length(q_range)
    @inbounds for i in 1:L
        q = q_range[i]
        ζ, S = zeta_estimator(X, q, τ_range)
        buffer[i, 1] = ζ
        buffer[i, 2] = S
    end
end

function zeta_estimator_range(X; τ_range = 1:19, q_range = 0.:0.1:1.)
    L = length(q_range)
    buffer = Matrix{Float64}(undef, L, 2)
    zeta_estimator_range!(buffer, X; τ_range = τ_range, q_range = q_range)
    buffer
end

function qth_abs_moment(X, τ, q)
    #check τ is feasible 
    if τ >= length(X)
        error("τ is too large!")
    end
    #Q = E[|X(t+ τ) - X(t)|^q]
    Q = 0.
    C = CartesianIndices(X)
    L = last(C)
    I1 = oneunit(L)
    #with gap τ
    for t in first(C):(last(C) - (τ * I1))
        @inbounds Q += abs(X[t + (τ * I1)] - X[t]) ^ q
    end
    Q / (length(C) - τ)
end

function generalised_hurst_range(X; τ_range = 1:19, q_range = 0.:0.1:1.)
    L = length(q_range)
    out = Matrix{Float64}(undef, L, 2)
    zeta_estimator_range!(out, X, τ_range = τ_range, q_range = q_range)
    @. out[:, 1] = out[:, 1] / q_range
    out
end

function generalised_hurst_range!(buffer, X; τ_range = 1:19, q_range = 0.:0.1:1.)
    L = length(q_range)
    zeta_estimator_range!(buffer, X, τ_range = τ_range, q_range = q_range)
    @. buffer[:, 1] = buffer[:, 1] / q_range
end

function generalised_hurst_exponent(X, q; τ_range = 1:19)
    generalised_hurst_range(X; τ_range = τ_range, q_range = q)
end

function hurst_exponent(X; τ_range = 1:19)
    generalised_hurst_range(X, τ_range = τ_range, q_range = 1)
end
