#=

    profile.jl

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me
=#

""" 
    profile(loglike,θ₀,lb,ub;kwargs...)

Profile the likelihood function `loglike`. θ₀ should generally be given as the (previously calculated) MLE for efficiency.

Required inputs `lb` and `ub` refer to the parameter bounds. Optional keyword arguments `lb_p` and `ub_p` can restrict the region that is profiled.

Example:

    loglike = θ -> loglikelihood(...,θ)
    
    Ψ,P  = profile(loglike,θ₀,lb,ub)
    plot(Ψ[1],P[1])

"""
function profile(loglike,θ₀,lb,ub;
        lb_p = lb,
        ub_p = ub,
        npts=100,
        params=1:length(θ₀)
    )

    dims = length(θ₀)

    Ψ = collect.(range.(lb_p,ub_p,length=npts))
    P = Array{Any,1}(undef,dims)

    # Loop through variables
    @threads for i = params

        P[i] = fill(-Inf,npts)

        lb1 = lb[setdiff(1:dims,i)]
        ub1 = ub[setdiff(1:dims,i)]

        # Start close to the MLE
        ψopt = θ₀[i]
        above_idx = Ψ[i][end] > ψopt ? findfirst(Ψ[i] .> ψopt) : length(Ψ[i]) + 1

        # Profile above the MLE
        θ01 = θ₀[setdiff(1:dims,i)]
        for j = above_idx:length(Ψ[i])
            ψ = Ψ[i][j]
            θ01,P[i][j] = optimise(λ -> loglike([ψ;λ][invperm([i;setdiff(1:dims,i)])]), θ01, lb1, ub1)
        end

        # Profile below the MLE
        θ01 = θ₀[setdiff(1:dims,i)]
        for j = above_idx-1:-1:1
            ψ = Ψ[i][j]
            θ01,P[i][j] = optimise(λ -> loglike([ψ;λ][invperm([i;setdiff(1:dims,i)])]), θ01, lb1, ub1)
        end

    end

    return Ψ,P

end