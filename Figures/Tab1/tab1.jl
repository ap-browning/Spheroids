#=
    Table 1

    Steady state parameter estimates and confidence intervals

    Runtime: Approximately 10 seconds (after module compiled)

=#

using CSV
using DataFrames
using DataFramesMeta
using Latexify

using Greenspan
using Inference

#########################################################################
## Load and subset data

    # Cell line and condition we are interested in
    CellLine            = "983b"
    InitialConditions   = [2500.0,5000.0,10000.0]

    # Look at day 18 (5000 and 10000) and day 21 (2500)
    confocal = @subset(CSV.read("Data/AllConfocalData.csv",DataFrame),
        :CellLine .== CellLine,
        (:Day .== 21.0) .& (:InitialCondition .== 2500) .|
        (:Day .== 18.0) .& (:InitialCondition .== 5000) .|
        (:Day .== 18.0) .& (:InitialCondition .== 10000))

    # Observation process with pooled covariance
    Σ = pooled_cov(confocal,:InitialCondition,[:R,:ϕ,:η])
    obs_process = M -> MvNormal(M,Σ)

#########################################################################
## Observation means: Compute MLEs, confidence intervals and perform tests

    # Model
    model = M -> M

    # Parameter bounds and initial guess
    lb    = [ 10.0, 0.0, 0.0]
    ub    = [500.0, 1.0, 1.0]
    θ₀    = [300.0, 0.8, 0.7]

    # Find MLEs and confidence intervals for each initial condition
    θmle1,CI1 = [Array{Any,1}(undef,3) for i = 1:3]
    θmle2,CI2 = [Array{Any,1}(undef,3) for i = 1:3]
    for (i,InitialCondition) ∈ enumerate(InitialConditions)

        # Log-likelihood function
        local loglike = θ -> loglikelihood(
            (obs_process ∘ model)(θ),   # Xᵢ ∼ f(m(θ))
            @subset(confocal, :InitialCondition .== InitialCondition)[:,[:R,:ϕ,:η]])

        # Find MLE and CI
        θmle1[i],CI1[i] = find_mle(loglike,θ₀,lb,ub)      

    end

    # Test for θ₂₅₀₀ = θ₅₀₀₀ 
    model = M -> M
    res11 = twosampletest(model,obs_process,θ₀,
        @subset(confocal, :InitialCondition .== 2500)[:,[:R,:ϕ,:η]],
        @subset(confocal, :InitialCondition .== 5000)[:,[:R,:ϕ,:η]];lb,ub
    )

    # Test for θ₅₀₀₀ = θ₁₀₀₀₀
    model = M -> M
    res12 = twosampletest(model,obs_process,θ₀,
        @subset(confocal, :InitialCondition .== 5000)[:,[:R,:ϕ,:η]],
        @subset(confocal, :InitialCondition .== 10000)[:,[:R,:ϕ,:η]];lb,ub
    )

#########################################################################
## Greenspan model: Compute MLEs, confidence intervals and perform tests

    # Model
    model = θ -> steady_state(θ)

    # Parameter bounds and initial guess
    lb    = [0.01,  10.0, 0.01]
    ub    = [0.99, 300.0, 50.00]
    θ₀    = [0.90, 250.0, 3.0]

    # Find MLEs and confidence intervals for each initial condition
    θmle2,CI2 = [Array{Any,1}(undef,3) for i = 1:3]
    for (i,InitialCondition) ∈ enumerate(InitialConditions)

        # Log-likelihood function
        local loglike = θ -> loglikelihood(
            (obs_process ∘ model)(θ),   # Xᵢ ∼ f(m(θ))
            @subset(confocal, :InitialCondition .== InitialCondition)[:,[:R,:ϕ,:η]])

        # Find MLE and CI
        θmle2[i],CI2[i] = find_mle(loglike,θ₀,lb,ub)

    end

    # Test for θ₂₅₀₀ = θ₅₀₀₀ 
    res21 = twosampletest(model,obs_process,θ₀,
        @subset(confocal, :InitialCondition .== 2500)[:,[:R,:ϕ,:η]],
        @subset(confocal, :InitialCondition .== 5000)[:,[:R,:ϕ,:η]];lb,ub
    )

    # Test for θ₅₀₀₀ = θ₁₀₀₀₀
    res22 = twosampletest(model,obs_process,θ₀,
        @subset(confocal, :InitialCondition .== 5000)[:,[:R,:ϕ,:η]],
        @subset(confocal, :InitialCondition .== 10000)[:,[:R,:ϕ,:η]];lb,ub
    )

#########################################################################
## Create table and convert to Latex
sigdigits = 3

tab1 = DataFrame(
    "Parameter"   => ["R", "ϕ", "η", "M", "Q", "γ", "Rcrit", "θ"],
    "Est (2500)"  => [round.(θmle1[1];sigdigits); NaN; round.(θmle2[1];sigdigits); NaN],
    "CI (2500)"   => [[round.(CI1[1][i];sigdigits) for i = 1:3]; NaN; [round.(CI2[1][i];sigdigits) for i = 1:3]; NaN],
    "Est (5000)"  => [round.(θmle1[2];sigdigits); NaN; round.(θmle2[2];sigdigits); NaN],
    "CI (5000)"   => [[round.(CI1[2][i];sigdigits) for i = 1:3]; NaN; [round.(CI2[2][i];sigdigits) for i = 1:3]; NaN],
    "Est (10000)" => [round.(θmle1[3];sigdigits); NaN; round.(θmle2[3];sigdigits); NaN],
    "CI (10000)"  => [[round.(CI1[3][i];sigdigits) for i = 1:3]; NaN; [round.(CI2[3][i];sigdigits) for i = 1:3]; NaN],
    "p (2500 to 5000)"  => [round.(res11.Pval;sigdigits=5); round.(res21.Pval;sigdigits=5)],
    "p (5000 to 10000)" => [round.(res12.Pval;sigdigits=5); round.(res22.Pval;sigdigits=5)]
)


res = latexify(tab1,env=:tabular)
write("$(@__DIR__)/tab1.tex",res)