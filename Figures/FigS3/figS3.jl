#=
    Figure 6

    Model profile likelihoods and 3D confidence regions at "steady state"

    Runtime: Approximately 40 seconds (after module compiled)

=#

using CSV
using DataFrames
using DataFramesMeta
using Distributions
using LinearAlgebra
using Statistics
using Plots
using Revise
using Greenspan
using Inference

include("../FigureDefaults.jl")

#########################################################################
## Load and subset data

    # Cell line and condition we are interested in
    CellLine            = "983b"
    InitialConditions   = [2500.0,5000.0,10000.0]

    # Look at day 18 (5000 and 10000) and day 21 (2500)
    confocal = @subset(CSV.read("Data/ConfocalData.csv",DataFrame),
        :CellLine .== "983b",
        (:Day .== 21.0) .& (:InitialCondition .== 2500) .|
        (:Day .== 21.0) .& (:InitialCondition .== 5000) .|
        (:Day .== 21.0) .& (:InitialCondition .== 10000))

    # Observation process with pooled covariance
    Σ = pooled_cov(confocal,:InitialCondition,[:R,:ϕ,:η])
    obs_process = M -> MvNormal(M,Σ)

########################################################s#################
## Profiles and confidence regions for Greenspan model

    # Model
    model = θ -> steady_state(θ)
    model_dv = θ -> steady_state_dv(θ) # Return Jacobian for confidence regions

    # Parameter bounds and initial guess
    lb    = [0.01,  10.0,  0.01]
    ub    = [0.999, 400.0, 50.00]
    θ₀    = [0.90, 250.0,   3.0]

    # Region to profile
    lb_p  = [0.65,  0.0, 0.4]
    ub_p  = [0.999, 200.0, 1.2]

    # Find MLE, profiles and confidence region for each initial density
    θmle,Lmle,Ψ,P,CR  = [Array{Any,1}(undef,3) for i = 1:5]

    # Runtime: Approximately 20 seconds
    for (i,InitialCondition) ∈ enumerate(InitialConditions)

        # Log-likelihood function
        local loglike = θ -> loglikelihood(
            (obs_process ∘ model)(θ),   # Xᵢ ∼ f(m(θ))
            @subset(confocal, :InitialCondition .== InitialCondition)[:,[:R,:ϕ,:η]])

        # Find MLE
        θmle[i],Lmle[i] = optimise(loglike,θ₀,lb,ub)

        # Calculate profiles
        Ψ[i],P[i] = profile(loglike,θ₀,lb,ub;lb_p,ub_p)

        # Log-likelihood derivative
        local loglike_dv = θ -> loglikelihood_dv(θ,model_dv,obs_process,
                @subset(confocal, :InitialCondition .== InitialCondition)[:,[:R,:ϕ,:η]])

        # Calculate confidence region
        CR[i] = confidenceregion3D(θmle[i],loglike_dv,lb,ub)

    end

########################################################s#################
## Confidence regions for the observation means

    # Derivative of model
    model    = M -> M
    model_dv = M -> (M,(diagm ∘ ones ∘ size)(M))

    # Parameter bounds and initial guess
    lb    = [ 10.0, 0.0, 0.0]
    ub    = [500.0, 1.0, 1.0]
    θ₀    = [300.0, 0.8, 0.7]

    # Find MLE, profiles and confidence region for each initial density
    θmle_identity,CR_identity = [Array{Any,1}(undef,3) for i = 1:2]

    # Runtime: Approximately 15 seconds
    for (i,InitialCondition) ∈ enumerate(InitialConditions)

        # Log-likelihood function
        local loglike = θ -> loglikelihood(
                (obs_process ∘ model)(θ),   # Xᵢ ∼ f(M)
                @subset(confocal, :InitialCondition .== InitialCondition)[:,[:R,:ϕ,:η]])

        # Find MLE
        θmle_identity[i],_ = optimise(loglike,θ₀,lb,ub)

        # Log-likelihood derivative
        local loglike_dv = θ -> loglikelihood_dv(θ,model_dv,obs_process,
                @subset(confocal, :InitialCondition .== InitialCondition)[:,[:R,:ϕ,:η]])

        # Calculate confidence region
        CR_identity[i] = confidenceregion3D(θmle_identity[i],loglike_dv,lb,ub)

    end

#########################################################################
## Produce figure

# Plot profiles
plts_profiles = [plot(xlabel=varname) for varname ∈ ["Q","Rcrit","Gamma"]]
for v = 1:3
    for (i,InitialCondition) ∈ enumerate(InitialConditions)
        plot!(plts_profiles[v],Ψ[i][v],max.(-10.0,P[i][v] .- Lmle[i]),
            c      = colors[CellLine][InitialCondition],
            frange = -4.0,
            falpha = 0.1,
            lw     = 0.75,
            ylim   = (-3.5,0.5),
            label  = Int(InitialCondition)
        )
    end
    hline!(plts_profiles[v],[-quantile(Chisq(1),0.95) / 2], c=:black, ls=:dashdot, lw=2, label="95% CI")
end

# Confidence region plots (with shadows on axes)
plts_CR = [plot([],[],[],xlim=(275.0,375.0),ylim=(0.85,0.95),zlim=(0.5,0.95),legend=:none)  # Observation means
           plot([],[],[],xlim=(0.6,1.0),ylim=(50.0,220.0),zlim=(0.2,1.3),legend=:none)]    # Greenspan model

for (j,cr) ∈ enumerate([CR_identity,CR]), (i,InitialCondition) ∈ enumerate(InitialConditions)

    # Plot shadows
    plot_cr3D_proj!(plts_CR[j],cr[i];c=colors[CellLine][InitialCondition],α=0.1,pα=0.5,lw=0.75)

    # Plot regions
    plot_cr3D!(plts_CR[j],cr[i];c=colors[CellLine][InitialCondition],lw=2)

    plot!(plts_CR[j],
        camera=(40,50),
        xlabel=["R [um]","Q [-]"][j],
        ylabel=["Phi [-]","Rc [um]"][j],
        zlabel=["Eta [-]","Gamma [-]"][j];
        defaults_3D...
    )

end

figS3 = plot(plts_profiles...,plts_CR...,layout=@layout([grid(1,3); grid(1,2){0.6h}]),legend=:none,size=(700,550))
savefig(figS3, "$(@__DIR__)/figS3.svg")