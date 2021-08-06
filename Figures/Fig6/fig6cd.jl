#=
    Figure 6cd (793b)

    Behaviour of spheroid structure through time

    Runtime: Approximately 2 seconds (after module compiled)
    
=#

using CSV
using DataFrames
using DataFramesMeta
using Plots
using Roots
using StatsPlots

using Greenspan
using Inference

include("../FigureDefaults.jl")

#########################################################################
## Load and subset data

    # Cell line and condition we are interested in
    CellLine            = "793b"
    InitialConditions   = [2500.0,5000.0,10000.0]

    # Load all data to plot
    confocal_all = @subset(CSV.read("Data/ConfocalData.csv",DataFrame),:CellLine .== CellLine)

#########################################################################
## Obtain MLE using steady state data

    # Model
    model = θ -> steady_state(θ)

    # Parameter bounds and initial guess
    lb    = [0.01,  10.0,  0.01]
    ub    = [0.99, 400.0, 50.00]
    θ₀    = [0.90, 250.0,   3.0]

    # Data subset to calibrate
    confocal_subset = @subset(CSV.read("Data/ConfocalData.csv",DataFrame),
        (:Day .== 24.0) .& (:InitialCondition .== 2500) .|
        (:Day .== 24.0) .& (:InitialCondition .== 5000) .|
        (:Day .== 24.0) .& (:InitialCondition .== 10000))

    # Observation process with pooled covariance
    Σ = pooled_cov(confocal_subset,:InitialCondition,[:R,:ϕ,:η])
    obs_process = M -> MvNormal(M,Σ)

    # Find MLE for each density
    θmle = Array{Any,1}(undef,3)

    # Runtime: Approximately 1 second
    for (i,InitialCondition) ∈ enumerate(InitialConditions)

        # Log-likelihood function
        local loglike = θ -> loglikelihood(
            (obs_process ∘ model)(θ),   # Xᵢ ∼ f(m(θ))
            @subset(confocal_subset, :InitialCondition .== InitialCondition)[:,[:R,:ϕ,:η]])

        # Find MLE
        θmle[i],_ = optimise(loglike,θ₀,lb,ub)

    end

#########################################################################
## Plot data and Greenspan model

fig6c = @df confocal_all scatter(:R,:ϕ,:η,group=:InitialCondition,c = [colors[CellLine][d] for d ∈ [2500,5000,10000]],
            ms=4,
            size=(500,500),
            legend=:none,
            xlabel="R",ylabel="phi",zlabel="eta";
            defaults_3D...)

    # Calculate and plot structure model for each MLE
    for (i,InitialCondition) ∈ enumerate(InitialConditions)

        Rmax = steady_state(θmle[i])[1]
        R    = collect(range(0.0,Rmax,length=400))
        X    = [R hcat(structure_model(θmle[i][1:2],R)...)']
        plot!(fig6c,X[:,1],X[:,2],X[:,3],c = [colors[CellLine][d] for d ∈ [2500,5000,10000]][i])
        
    end

#########################################################################
## Subset data to perform ANOVA and test for dependence on initial condition

    # Subset
    confocal_all_1 = @subset(confocal_all,:η .> 0.2, :ϕ .> 0.68, :η .< 0.005 * :R .- 0.85)
    confocal_all_2 = @subset(confocal_all,(!).((:η .> 0.2) .& (:ϕ .> 0.68) .& (:η .< 0.005 * :R .- 0.85)))

    # "Orthogonal" linear model test
    p,β₀,β₁ = @df confocal_all_1 lm_orthogonal_test(:R,:ϕ,:η,:InitialCondition)
    bestfit = lm_orthogonal_linefcn(β₀)
    print("Likelihood ratio test p-value for $CellLine line: $p")

    # Plot best fit line (for 0 < η < 1)
    x₀ = bestfit(0.0)   
    x₁ = bestfit(find_zero(t -> bestfit(t)[3] - 1.0, 0.0)) # η = 1
    fig6d = plot([[x₀[i],x₁[i]] for i = 1:3]..., c=:black)

    @df confocal_all_1 scatter!(fig6d,:R,:ϕ,:η,group=:InitialCondition,c = [colors[CellLine][d] for d ∈ [2500,5000,10000]],
            ms=4,
            size=(500,500),
            legend=:none,
            xlabel="R",ylabel="phi",zlabel="eta";
            defaults_3D...)
    @df confocal_all_2 scatter!(fig6d,:R,:ϕ,:η,c=:black,α=0.2,ms=4)

#########################################################################
## Save figure 7

fig6cd = plot(fig6c, fig6d, link=:all, size=(800,400),
        xlim=(-10.0,410.0), 
        ylim=(-0.02,1.02),
        zlim=(-0.02,1.02) 
)

#savefig("$(@__DIR__)/fig6cd.svg")

