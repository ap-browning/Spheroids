#=
    Figure 3

    Confidence regions and profile likelihood example.

    Runtime: Approximately 20 seconds (after module compiled)
    
=#

using CSV
using DataFrames
using DataFramesMeta
using Distributions
using Plots
using Statistics
using StatsPlots

using Greenspan
using Inference

include("../FigureDefaults.jl")

#########################################################################
## Load and subset data

    # Look at 983b, 5000, D14
    confocal = @subset(CSV.read("Data/ConfocalData.csv",DataFrame),
                    :CellLine .== "983b",
                    :InitialCondition .== 5000,
                    :Day .== 14.0)

#########################################################################
## Model and likelihood

    # Get mean outer radius
    R = @df confocal mean(:R)

    # Get standard deviation for observation process
    Σ = cov(Array(confocal[:,[:ϕ,:η]]))
    obs_process = M -> MvNormal(M,Σ)

    # Model
    model = ξ -> structure_model(ξ,R)
    model_dv = ξ -> structure_model_dv(ξ,R)

    # Likelihood function
    loglike = θ -> loglikelihood((obs_process ∘ model)(θ),confocal[:,[:ϕ,:η]])
    loglike_dv = θ -> loglikelihood_dv(θ,model_dv,obs_process,confocal[:,[:ϕ,:η]])

    # Parameter bounds
    lb = [0.10, 10.0]
    ub = [0.99,500.0]
    θ₀ = [0.80,200.0]

    # Get MLE
    θmle,Lmle = optimise(loglike,θ₀,lb,ub)

    # Get 50% confidence region
    cr50 = confidenceregion2D(θmle,loglike_dv;α=0.5)
    cr95 = confidenceregion2D(θmle,loglike_dv;α=0.05)

#########################################################################
## Plot likelihood surface

    Qrange = range(0.655,0.735,length=300)
    Rrange = range(184.0,201.0,length=300)

    Qticks = 0.66:0.01:0.73
    Rticks = 185:5:200

    L = [loglike([Q,R₁]) for Q ∈ Qrange, R₁ ∈ Rrange]
    L .-= maximum(L)  

    plt = plot(Qrange,Rrange,L',st=:contourf,lw=0.0,clim=(-8.0,0.0),c=:RdPu_9,extendlow=-8)

    plot!(plt,cr50[:,1],cr50[:,2],lw=3.0,c=:white,ls=:dot)
    plot!(plt,cr95[:,1],cr95[:,2],lw=3.0,c=:white,ls=:dash)

    plot!(plt,
        axis=:all, box=:on,legend=:none,widen=:false,grid=:off,
        xticks=Qticks,yticks=Rticks,xlabel="Q [um]",ylabel="Rc [um]"
    )

#########################################################################
## Profiles

    # Profiles (manually to get trace)
    Λ1 = zeros(length(Qrange))
    L1 = zeros(length(Qrange))
    for i = 1:length(Qrange)
        λ,L1[i] = optimise(R₁ -> loglike([Qrange[i];R₁]), [θmle[2]], [lb[2]], [ub[2]])
        Λ1[i] = λ[1]
    end

    Λ2 = zeros(length(Rrange)) 
    L2 = zeros(length(Rrange))
    for i = 1:length(Rrange)
        λ,L2[i] = optimise(Q -> loglike([Q; Rrange[i]]), [θmle[1]], [lb[1]], [ub[1]])
        Λ2[i] = λ[1]
    end

    # Add profiles to 2D plot
    plot!(plt,Qrange,Λ1,c=:black)
    plot!(plt,Λ2,Rrange,c=:black,ls=:dash)

    # Plot MLE
    scatter!(plt,[θmle[1]],[θmle[2]],c=:white)

    # Plot profile 1
    inconfint_1 = L1 .- maximum(L1) .> -quantile(Chisq(1),0.95)/2
    plt_Q = plot(Qrange,L1 .- maximum(L1),c=:black,α=0.2,legend=:none,ylim=(-3.1,0.1),xlabel="Q [um]",xticks=Qticks,widen=:false)
    plot!(plt_Q,Qrange[inconfint_1],L1[inconfint_1] .- maximum(L1),c=:black)
    hline!(plt_Q,-[quantile(Chisq(1),0.95)/2],c=:red)

    # Plot profile 2
    inconfint_2 = L2 .- maximum(L2) .> -quantile(Chisq(1),0.95)/2
    plt_R = plot(Rrange,L2 .- maximum(L2),c=:black,α=0.2,legend=:none,ylim=(-3.1,0.1),xlabel="Rc [um]",ls=:dash,xticks=Rticks,widen=:false)
    plot!(plt_R,Rrange[inconfint_2],L2[inconfint_2] .- maximum(L2),c=:black,ls=:dash)
    hline!(plt_R,-[quantile(Chisq(1),0.95)/2],c=:red)

#########################################################################
## Plot

fig3 = plot(plt,plt_Q,plt_R,layout=@layout([a{0.6h}; grid(2,1)]),size=(600,800))
#savefig(fig3,"$(@__DIR__)/fig3.svg")
