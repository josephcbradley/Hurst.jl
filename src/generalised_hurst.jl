import LinearAlgebra: dot

function qth_abs_moment(X::Vector{T}, τ, q) where {T<:Real}
    #Calculte Q = E[|X(t+ τ) - X(t)|^q]
    #check τ is feasible 
    if τ <= 0 || τ >= length(X)
        throw(DimensionMismatch("τ is out of bounds!"))
    end
    #Q = E[|X(t+ τ) - X(t)|^q]
    Q = 0.0
    n = length(X)
    #with gap τ
    for t = 1:(n-τ)
        @inbounds Q += abs(X[t+τ] - X[t])^q
    end
    Q / (n - τ)
end

function lsq_estimator(Y::Vector{T}, S::Vector{T}) where {T<:Real}
    N = length(Y)
    if N != length(S)
        throw(DimensionMismatch("Dimensions do not match!"))
    end
    Ss = sum(S)
    Sy = sum(Y)
    Ssy = dot(S, Y)
    Sss = dot(S, S)
    Syy = dot(Y, Y)

    # Check for division by zero
    den = N * Sss - Ss^2
    if N <= 2 || den <= 0
        return (NaN, NaN)
    end

    ζ = (N * Ssy - Ss * Sy) / den
    sϵ2 = (1 / (N * (N - 2))) * (N * Syy - Sy^2 - ζ^2 * den)
    SD = sqrt((N * sϵ2) / den)
    ζ, SD
end

function zeta_estimator!(
    Y::Vector{T},
    S::Vector{T},
    X::Vector{T},
    τ_range,
    q;
    use_threading::Bool=true,
) where {T<:Real}
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

    if (length(Y) != N) || (length(S) != N)
        throw(DimensionMismatch("Dimensions do not match!"))
    end

    #calculate regression data
    if use_threading
        Threads.@threads for i = 1:N
            @inbounds begin
                τ = τ_range[i]
                Y[i] = log(qth_abs_moment(X, τ, q))
                S[i] = log(τ)
            end
        end
    else
        for i = 1:N
            @inbounds begin
                τ = τ_range[i]
                Y[i] = log(qth_abs_moment(X, τ, q))
                S[i] = log(τ)
            end
        end
    end

    #simple regression formulae
    ζ, SD = lsq_estimator(Y, S)

    #note ζ(q) = qH(q), so H(q) = ζ(q) / q 
    return ζ, SD
end

function zeta_estimator(X, τ_range, q; use_threading::Bool=true)
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
    zeta_estimator!(Y, S, X, τ_range, q; use_threading=use_threading)
end

function zeta_estimator_range!(range_buffer, Y, S, X, τ_range, q_range; use_threading::Bool=true)
    L = length(q_range)
    if size(range_buffer) != (L, 2)
        throw(DimensionMismatch("range_buffer has incorrect dimensions"))
    end
    @inbounds for i = 1:L
        q = q_range[i]
        ζ, SD = zeta_estimator!(Y, S, X, τ_range, q; use_threading=use_threading)
        range_buffer[i, 1] = ζ
        range_buffer[i, 2] = SD
    end
end

"""
    zeta_estimator_range(X, τ_range, q_range)

Calculate ``\\zeta (q)`` that satifies:

``\\ln \\left(E\\left[|X(t+\\tau)-X(t)|^q\\right]\right)=\\zeta(q) \\ln (\\tau)+\\ln (K(q))``

for some series ``X(t)``, over the vector `q_range`. 

Returns a `(length(q_range), 2)` matrix where the first column contains the values of the ``\\zeta(q)`` for different `q` and the second column contains the standard errors.

See also [`hurst_exponent`](@ref).
"""
function zeta_estimator_range(X, τ_range, q_range; use_threading::Bool=true)
    L = length(q_range)
    range_buffer = Matrix{Float64}(undef, L, 2)
    N = length(τ_range)
    Y = Vector{Float64}(undef, N)
    S = similar(Y)
    zeta_estimator_range!(range_buffer, Y, S, X, τ_range, q_range; use_threading=use_threading)
    range_buffer
end

function generalised_hurst_range!(buffer, Y, S, X, τ_range, q_range; use_threading::Bool=true)
    L = length(q_range)
    if size(buffer) != (L, 2)
        throw(DimensionMismatch("buffer has incorrect dimensions"))
    end
    zeta_estimator_range!(buffer, Y, S, X, τ_range, q_range; use_threading=use_threading)
    # If the process does not actually scale in the assumed way, the linear regression method will throw values 
    # for H outside [0, 1]. To make this clear, we will set any such values of the buffer to NaN 
    for i = 1:L
        @inbounds begin
            q = q_range[i]
            H = buffer[i, 1] / q
            if isfinite(H) && 0 <= H <= 1
                buffer[i, 1] = H
            else
                buffer[i, 1] = NaN
                buffer[i, 2] = NaN
            end
        end
    end
end

"""
    generalised_hurst_range(X, τ_range, q_range)

Calculate the generalised Hurst exponent (GHE) of the series `X` with absolute moments `q_range` over the range `τ_range`, along with its standard error.

Returns a `(length(q_range), 2)` matrix where the first column contains the values of the GHE and the second column contains the standard errors.

See also [`hurst_exponent`](@ref).

# Examples
```jldoctest
julia> X = accumulate(+, randn(1000));

julia> q_range = 0.:0.1:1.; tau_range = 1:19;

julia> generalised_hurst_range(X, tau_range, q_range);
```
"""
function generalised_hurst_range(X, τ_range, q_range; use_threading::Bool=true)
    L = length(q_range)
    out = Matrix{Float64}(undef, L, 2)
    N = length(τ_range)
    Y = Vector{Float64}(undef, N)
    S = similar(Y)

    generalised_hurst_range!(out, Y, S, X, τ_range, q_range; use_threading=use_threading)
    out
end

"""
    generalised_hurst_exponent(X, τ_range, q)

Calculate the generalised Hurst exponent of the series `X` with absolute moment `q` over the range `τ_range` along with its standard error.

See also [`hurst_exponent`](@ref).

# Examples
```jldoctest
julia> X = accumulate(+, randn(1000));

julia> generalised_hurst_exponent(X, 1:19, 0.5);
```
"""
generalised_hurst_exponent(X, τ_range, q; use_threading::Bool=true) = generalised_hurst_range(X, τ_range, q; use_threading=use_threading)

"""
    hurst_exponent(X, τ_range)

Calculate the Hurst exponent of the series `X` over the range `τ_range` along with its standard error.

See [Buonocore et al. 2016](https://doi.org/10.1016/j.chaos.2015.11.022).

# Examples
```jldoctest
julia> X = accumulate(+, randn(1000));

julia> isapprox(hurst_exponent(X, 1:19)[1], 0.5, atol = 0.1)
true
```
"""
function hurst_exponent(X::Vector{T}, τ_range; use_threading::Bool=true) where {T<:Real}
    generalised_hurst_range(X, τ_range, 1; use_threading=use_threading)
end
