module Hurst

import GLM, DataFrames

# Write your package code here.
include("generalised_hurst.jl")

export 
    qth_abs_moment,
    zeta_estimator,
    generalised_hurst_range, 
    generalised_hurst_exponent, 
    hurst_exponent, 
    zeta_estimator_range
    generalised_hurst_range!
end
