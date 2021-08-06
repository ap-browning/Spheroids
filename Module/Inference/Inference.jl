#=

    Inference.jl

    A Julia module to perform inference, including profile likelihood, confidence 
    regions, etc.

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=#
module Inference

    using DataFrames
    using DataFramesMeta
    using DifferentialEquations
    using Distributions
    using DiffResults
    using ForwardDiff
    using Interpolations
    using LazySets
    using LinearAlgebra
    using NLopt
    using Plots
    using Roots
    using StatsBase
    using StatsFuns
    using .Threads

    import StatsBase: loglikelihood

    export optimise, profile, find_mle, twosampletest, loglikelihood, loglikelihood_dv,
        confidenceregion2D, confidenceregion3D, lm_orthogonal, lm_orthogonal_linefcn, lm_orthogonal_test,
        plot_cr3D, plot_cr3D!, plot_cr3D_proj!, pooled_data, pooled_cor, pooled_cov, pooled_std

    include("confidence_regions.jl")
    include("maxlikelihood.jl")
    include("optimise.jl")
    include("pool_data.jl")
    include("profile.jl")
    include("total_least_squares.jl")
    
end