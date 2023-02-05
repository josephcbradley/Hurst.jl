module Hurst

import LinearAlgebra

# Write your package code here.
include("generalised_hurst.jl")

export hurst_exponent, generalised_hurst_exponent, generalised_hurst_range
zeta_estimator_range
end
