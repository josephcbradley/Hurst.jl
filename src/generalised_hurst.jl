function hurst_estimator(X, τ_range = 1:19)
    #ζ(q) = qH(q) where H is the GHE 
    #let Q = E[|X(t+ τ) - X(t)|^q] = K(q)τ^ζ(q)
    #take logarithms and rearrange
    #ln(Q) = ζ(q)ln(τ) + ln(K(q))
    #to estimate τ, consider the linear regression of
    #Y = β₀ + β₁S where 
    #Y = ln(Q), S = ln(τ), β₀ = ln(K(q)), and β₁ = ζ(q)
    #return β₁ and its std error

    #form Y 
    N = length(τ_range)
    Y = Vector{Float64}(undef, N)
    β₀ = Vector{Float64}(undef, N)
    β₁ = Vector{Float64}(undef, N)
    S = Vector{Float64}(undef, N)

    for i in CartesianIndices(τ_range)
        τ = τ_range[i]
        Q = qth_abs_moment(X, τ)
    end

end

function qth_abs_moment(X, τ, q)
    out = 0.
    C = CartesianIndices(X)
    L = last(C)
    I1 = oneunit(L)
    #with gap τ, start counting after 1 + τ indices
    for t in (1 + I1*τ):L
        out += (abs(X[t] - X[t - I1*τ]) / abs(X[t - I1*τ])) ^ q
    end
    out
end