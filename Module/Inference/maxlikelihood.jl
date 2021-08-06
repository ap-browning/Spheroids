#=

    maxlikelihood.jl

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=#

loglikelihood(d,x::DataFrame) = loglikelihood(d,Array(x)')

""" 
    loglikelihood_dv(θ,model_dv,obs_process,data)

Computes the loglikelihood and the loglikelihood gradient.

Example:

    using Greenspan
    model_dv = θ -> steady_state_dv(θ)
    obs_process = M -> MvNormal(M,Σ)

    loglike_dv = θ -> loglikelihood_dv(θ,model_dv,obs_process,data)

"""
function loglikelihood_dv(θ::Vector,model_dv::Function,obs_process::Function,data)
    M,J = model_dv(θ)
    result = DiffResults.GradientResult(M)
    ForwardDiff.gradient!(result, M -> loglikelihood(obs_process(M),data), M)
    DiffResults.value(result), J' * DiffResults.gradient(result)
end


""" 
    find_mle(loglike,θ₀,lb,ub;kwargs...)

Finds maximum likelihood estimates and confidence intervals directly using profile likelihood.

Example:

    loglike = θ -> loglikelihood(...,θ)
    
    mle,ci95 = find_mle(loglike,θ₀,lb,ub)

"""
function find_mle(loglike,θ₀,lb,ub;α=0.95,xrtol=1e-3)

    dimθ = length(θ₀)

    ## Calculate maximum likelihood estimate
    θmle,Lmle = optimise(loglike,θ₀,lb,ub)

    ## Loop through parameters and find confidence intervals using "profile likelihood" method

        CI = Array{Array{Float64,1},1}(undef,dimθ)
        
        # Log likelihood threshold
        thresh = Lmle - quantile(Chisq(1),α) / 2

        # Loop through parameters
        @threads for idx_ψ = 1:dimθ

            # Nuisance parameter indices
            idx_λ = setdiff(1:dimθ,idx_ψ)

            # Sectioned likelihood
            loglike1 = (ψ,λ) -> loglike([ψ;λ][invperm([idx_ψ;idx_λ])])

            # Find upper and lower values of ψ such that max_λ loglike1(ψ) > thresh
            maxloglike1 = ψ -> optimise(λ -> loglike1(ψ,λ),θmle[idx_λ],lb[idx_λ],ub[idx_λ])[2]
            lower = maxloglike1(lb[idx_ψ]) > thresh ? lb[idx_ψ] : find_zero(ψ -> maxloglike1(ψ) - thresh,(lb[idx_ψ],θmle[idx_ψ]), Roots.A42())
            upper = maxloglike1(ub[idx_ψ]) > thresh ? ub[idx_ψ] : find_zero(ψ -> maxloglike1(ψ) - thresh,(θmle[idx_ψ],ub[idx_ψ]), Roots.A42())

            CI[idx_ψ] = [lower,upper]

        end

    θmle,CI

end

"""
    twosampletest(model,obs_process,x₁,x₂,θ₀,...)

Conduct a two sample hypothesis test for parameter equivalence.

Example:

    model       = θ -> steady_state(θ)
    obs_process = M -> MvNormal(M,Σ)
    θ₀          = rand(Product(Uniform.(lb,ub)))
        
    testres     = twosampletest(model,obs_process,data₁,data₂,θ₀)

Note:
    `x₁` and `x₂` may be a data frame (with number columns equal to `length(M)`) or a `Matrix`.
    
"""
function twosampletest(model::Function,obs_process::Function,θ₀,x₁,x₂;
        lb      = fill(0.0,size(θ₀)),
        ub      = fill(Inf,size(θ₀)),
        names   = ["p$i" for i = 1:length(θ₀)],
        α       = 0.05
    )

    dimθ = length(θ₀)

    ## Test whether θ̂₁ = θ̂₂

        # Calculate log-likelihood values (θ̂₁ ≂̸ θ̂₂)
        θmle₁,Lmle₁ = optimise(θ -> loglikelihood((obs_process ∘ model)(θ),x₁),θ₀,lb,ub)
        θmle₂,Lmle₂ = optimise(θ -> loglikelihood((obs_process ∘ model)(θ),x₂),θ₀,lb,ub)

        # Calculate log-likelihood values (θ̂₁ = θ̂₂)
        _,Lmle₀ = optimise(θ -> loglikelihood((obs_process ∘ model)(θ),[x₁ x₂]),θ₀,lb,ub)

    # Test statistic, λ ∼ χ²(dim(θ))
    λ     = 2(Lmle₁ + Lmle₂ - Lmle₀)
    p     = 1.0 - cdf(Chisq(dimθ),λ)

    ## Test whether θ₁ⁱ = θ₂ⁱ
    function funp(θp,i)

        # Unpack θp -> θ₁ and θ₂ (i.e., models are only constrained to have ith parameters equal)
        # θp = [θⁱ,θ̃₁,θ̃₂]
        θ₁,θ₂ = zeros(dimθ),zeros(dimθ)
        θ₁[i] = θp[1]
        θ₂[i] = θp[1]
        θ₁[setdiff(1:dimθ,i)] = θp[2:2+dimθ-2]
        θ₂[setdiff(1:dimθ,i)] = θp[2+dimθ-1:end]

        loglikelihood((obs_process ∘ model)(θ₁),x₁) + loglikelihood((obs_process ∘ model)(θ₂),x₂)

    end
    pvals = zeros(dimθ)
    λvals = zeros(dimθ)
    for i = 1:dimθ

        ni = setdiff(1:dimθ,i)
        _,Lmleᵢ = optimise(θ -> funp(θ,i),[θ₀[i]; θ₀[ni]; θ₀[ni]],[lb[i]; lb[ni]; lb[ni]],[ub[i]; ub[ni]; ub[ni]])

        λvals[i] = 2(Lmle₁ + Lmle₂  - Lmleᵢ)
        pvals[i] = 1.0 - cdf(Chisq(1),λvals[i])

    end

    DataFrame(
        Test    = [names .* "₁ = " .* names .* "₂";"θ₁ = θ₂"],
        Est1    = [θmle₁; missing],
        Est2    = [θmle₂; missing],
        Score   = [λvals; λ],
        df      = [ones(Int,dimθ); dimθ],
        Pval    = [pvals; p]
    )

end
function twosampletest(model::Function,obs_process::Function,θ₀,x₁::DataFrame,x₂::DataFrame;kwargs...)
    twosampletest(model,obs_process,θ₀,Array(x₁)',Array(x₂)';kwargs...)
end