#=

    Greenspan.jl

    A Julia module to solve the Greenspan ODE model and the associate structure and steady state problems.

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me


                Jesse A. Sharp
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                jesse.sharp@hdr.qut.edu.au

=#
module Greenspan

    using DifferentialEquations
    using DiffResults
    using ForwardDiff
    using Polynomials
    using Roots

    export coefs, steady_state, steady_state_dv, solve_transient_model, 
        solve_transient_model_all_vars, solve_transient_model_all_vars_dim,
        structure_model, structure_model_dv

    include("steady_state.jl")
    include("structure_model.jl")
    include("transient_model.jl")
    
end