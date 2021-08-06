#=
    Figure 4

    Observation means and parameter estimates using data at each time point.

    Runtime: Approximately 20 seconds (after module compiled)
    
=#

using CSV
using DataFrames
using DataFramesMeta
using Distributions
using Statistics
using Plots

using Greenspan
using Inference

include("../FigureDefaults.jl")

#########################################################################
## Load and subset data

    # Cell line and condition we are interested in
    CellLine            = "983b"
    InitialConditions   = [2500.0,5000.0,10000.0]

    # Look at day 12 data onwards
    confocal = @subset(CSV.read("Data/ConfocalData.csv",DataFrame),:CellLine .== CellLine, :Day .≥ 12.0)

#########################################################################
## "Identity model": fit the mean of each observation at each time

    # Parameter bounds and initial guess
    lb    = [ 10.0, 0.0, 0.0]
    ub    = [500.0, 1.0, 1.0]
    θ₀    = [300.0, 0.8, 0.7]

    # Find MLE and 95% confidence interval at each day
    θmle1 = Array{Any,1}(undef,length(InitialConditions))
    CI1   = Array{Any,1}(undef,length(InitialConditions))
    Days1 = Array{Any,1}(undef,length(InitialConditions))

    # Runtime: Approximately 10 seconds
    for (i,InitialCondition) ∈ enumerate(InitialConditions)
        Days1[i] = sort(unique(@subset(confocal,:InitialCondition .== InitialCondition).Day))
        θmle1[i] = Array{Any,1}(undef,length(Days1[i]))
        CI1[i]   = Array{Any,1}(undef,length(Days1[i]))

        for (j,Day) ∈ enumerate(Days1[i])

            # Subset of data
            local confocal_subset = @subset(confocal,:InitialCondition .== InitialCondition, :Day .== Day)

            # Covariance matrix
            local Σ = cov(Array(confocal_subset[:,[:R,:ϕ,:η]]))

            # Log-likelihood
            local loglike = θ -> loglikelihood(MvNormal(θ,Σ),confocal_subset[:,[:R,:ϕ,:η]])

            # Find MLE and confidence intervals
            θmle1[i][j],CI1[i][j] = find_mle(loglike,θ₀,lb,ub)

        end

    end

#########################################################################
## Greenspan model: fit the steady state model at all time points. 
#  This is mathematically equivalent to fitting the structure model, and inferring the mean 
#  outer radius at each time point (for parameter Q and R_crit)

    model = θ -> steady_state(θ)

    # Parameter bounds and initial guess
    lb    = [0.10,  10.0,  0.01]
    ub    = [0.99, 400.0, 50.00]
    θ₀    = [0.90, 250.0,   3.0]

    # Find MLE and 95% confidence interval at each day
    θmle2 = Array{Any,1}(undef,length(InitialConditions))
    CI2   = Array{Any,1}(undef,length(InitialConditions))
    Days2 = Array{Any,1}(undef,length(InitialConditions))

    # Runtime: Approximately 10 seconds
    for (i,InitialCondition) ∈ enumerate(InitialConditions)
        Days2[i] = sort(unique(@subset(confocal,:InitialCondition .== InitialCondition).Day))
        θmle2[i] = Array{Any,1}(undef,length(Days2[i]))
        CI2[i]   = Array{Any,1}(undef,length(Days2[i]))

        for (j,Day) ∈ enumerate(Days2[i])

            # Subset of data
            local confocal_subset = @subset(confocal,:InitialCondition .== InitialCondition, :Day .== Day)

            # Covariance matrix
            local Σ = cov(Array(confocal_subset[:,[:R,:ϕ,:η]]))

            # Log-likelihood
            local loglike = θ -> loglikelihood(MvNormal(model(θ),Σ),confocal_subset[:,[:R,:ϕ,:η]])

            # Find MLE ansd confidence intervals
            θmle2[i][j],CI2[i][j] = find_mle(loglike,θ₀,lb,ub)

        end

    end

#########################################################################
## Produce figure

plts = [plot() for i = 1:2, j = 1:3]
θmle = [θmle1,θmle2]
CI   = [CI1,CI2]

for r = 1:2, c = 1:3, (i,InitialCondition) ∈ enumerate(InitialConditions)

    θ = [θmle[r][i][d][c]  for d = 1:length(Days1[c])]
    l = [CI[r][i][d][c][1] for d = 1:length(Days1[c])]
    u = [CI[r][i][d][c][2] for d = 1:length(Days1[c])]

    plot!(plts[r,c], Days1[i] .+ (i-2) * 0.2, θ, c = colors[CellLine][InitialCondition],
            lα=0.3, lw=1.0, ms=3, shape=:circle)
    for d = 1:length(Days1[c])
        plot!(plts[r,c],fill(Days1[c][d],2) .+ (i-2) * 0.2, [l[d],u[d]], c = colors[CellLine][InitialCondition],
            lw=1.5)
    end

    plot!(plts[r,c],
        title   = [["R [um]","Phi [-]","Eta [-]"],["Q [-]","Rc [um]","Gamma [-]"]][r][c],
        xticks  = [12.0,14.0,16.0,18.0,21.0],
        xlabel  = "Time [d]",
        ylim    = [[[270,380],[0.8,0.96],[0.4,0.9]],[[0.6,1.0],[70,250],[0,5]]][r][c]
    )

end

fig4 = plot(plts[1,:]...,plts[2,:]...,legend=:none,size=(800,450))
#savefig(fig4,"$(@__DIR__)/fig4.svg")
