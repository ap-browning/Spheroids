#=

    confidence_regions.jl

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me
=#

"""
    confidenceregion2D(ξopt,loglike,∇loglike)

Construct a bivariate confidence region for the likelihood function `loglike` around the MLE `ξopt`.

The gradient of the loglikelihood function must be given `∇loglike`. Alternatively, a single input `loglike_dv` can be given that returns both the likelihood value and gradient.

Example:

    model_dv    = θ -> structure_model_dv(θ)
    obs_process = M -> MvNormal(M,Σ)
    loglike_dv  = θ -> loglikelihood_dv(θ,model_dv,obs_process,data)
    
    θopt = optimise(θ -> loglikedv(θ)[1],θ₀,lb,ub)

    cr = confidenceregion2D(θopt,loglike_dv)
    plot(cr[:,1],cr[:,2])

"""
function confidenceregion2D(ξopt,loglike,∇loglike;
        α  = 0.05,
        lb = fill(0.0,length(ξopt)),
        ub = fill(Inf,length(ξopt)),
        method = Heun(),
        dt     = 0.1,
        tmax   = 1e3,
        τ      = loglike(ξopt) - quantile(Chisq(2),1.0 - α) / 2
    )

    # Bounds
    inbounds = u -> all(u .> lb) && all(u .< ub)

    # Start at the "top"
    d  = find_zero(d -> loglike(ξopt + [0.0,d]) - τ,(0.0, ub[2] - ξopt[2]))
    X₀ = ξopt + [0.0,d]

    # Move around the space perpendicular to ∇loglike
    function rhs(u,p,t)
        if inbounds(u)
            grad = ∇loglike(u)
            return [-grad[2] ; grad[1]] / norm(grad)
        else
            return [0.0; 0.0]
        end
    end

    # Stop when we've looped back
    isclockwise = ∇loglike(X₀)[2] < 0
    function callback(u,t,int)
        res = isclockwise ? (u[1] > X₀[1] > int.uprev[1]) : (u[1] < X₀[1] < int.uprev[1])
        return res
    end
    cb  = DiscreteCallback(
        callback,
        integrator -> terminate!(integrator))

    # Pause if we go out of bounds (todo: allow multiple out of bounds instances)
    cb2 = DiscreteCallback(
        (u,t,int) -> inbounds(u) == false,
        integrator -> terminate!(integrator))

    # Solve
    sol = solve(ODEProblem(rhs,X₀,(0.0,tmax)),method,dt=dt,callback=CallbackSet(cb,cb2),reltol=1e-5)
    X   = hcat(sol.u...)'

    # Did we go out of bounds?
    if inbounds(X[end,:]) == false

        display("Out of bounds once!")

        # Which bounds were violated
        violated = [X[end,:] .≤ lb X[end,:] .≥ ub]
        param    = findfirst(any(violated,dims=2))[2]

        # Last point within the bounds
        idx      = findlast([inbounds(X[i,:]) for i = 1:size(X,1)])
        X        = X[1:idx,:]
        Xlast    = X[end,:]
        Llast    = loglike(Xlast)

        # Which direction should we move in?
        gradlike = ∇loglike(Xlast)
        if param == 1 # Parameter 1 bounds, move parameter 2

            direction = sign(gradlike[2])
            if Llast < 0
                Xlast[2] = find_zero(p₂ -> loglike([Xlast[1],p₂]) - τ,Xlast[2])
            end
            lookat    = direction > 0 ? (Xlast[2]+1e-5,ub[2]) : (lb[2],Xlast[2]-1e-5)
            p₂        = find_zero(p₂ -> loglike([Xlast[1],p₂]) - τ,lookat)
            X₁        = [Xlast[1],p₂]

        else

            direction = sign(gradlike[1])
            if Llast < 0
                Xlast[1] = find_zero(p₁ -> loglike([p₁,Xlast[2]]) - τ,Xlast[1])
            end
            lookat    = direction > 0 ? (Xlast[1]+1e-5,ub[1]) : (lb[1],Xlast[1]-1e-5)
            p₁        = find_zero(p₁ -> loglike([p₁,Xlast[2]]) - τ,lookat)
            X₁        = [p₁,Xlast[2]]

        end

        sol = solve(ODEProblem(rhs,X₁,(0.0,tmax)),method,dt=dt,callback=CallbackSet(cb,cb2))
        X   = [X; hcat(sol.u...)']

    end

    X

end
function confidenceregion2D(ξopt,loglike_dv;kwargs...)
    confidenceregion2D(ξopt,θ -> loglike_dv(θ)[1], θ -> loglike_dv(θ)[2];kwargs...)
end


"""
    confidenceregion3D(ξopt,loglike,∇loglike)

Construct a trivariate confidence region for the likelihood function `loglike` around the MLE `θopt`.

The confidence region is constructed by taking slices (default: `slice = 20`) through the third (default: `slicevar = 3`) parameter.

The gradient of the loglikelihood function must be given `∇loglike`. Alternatively, a single input `loglike_dv` can be given that returns both the likelihood value and gradient.

Example:

    model_dv    = θ -> steady_state_dv(θ)
    obs_process = M -> MvNormal(M,Σ)
    loglike_dv  = θ -> loglikelihood_dv(θ,model_dv,obs_process,data)
    
    θopt = optimise(θ -> loglikedv(θ)[1],θ₀,lb,ub)

    cr3D = confidenceregion3D(θopt,loglike_dv)
    plotCR3D(cr3D)

"""
function confidenceregion3D(θopt,loglike,∇loglike,lb,ub;
        α      = 0.05,
        slice  = 20,
        method = Heun(),
        dt     = 0.1,
        tmax   = 1e3,
        slicevar = 3
    )

    # Likelihood value to contour
    τ  = loglike(θopt) - quantile(Chisq(3),1.0 - α) / 2

    # Variable permutationss
    idx_ξ = setdiff(1:3,slicevar)
    idx_λ = slicevar

    # Sectioned log likelihood
    loglike1 = (ξ,λ) -> loglike([ξ;λ][invperm([idx_ξ;idx_λ])])

    # Slices on auto or not?
    if typeof(slice) != Int
        P₃ = typeof(slice) == Float64 ? [slice] : slice
    else
        # Find upper and lower values of λ
        maxloglike1 = λ -> optimise(ξ -> loglike1(ξ,λ),θopt[idx_ξ],lb[idx_ξ],ub[idx_ξ])[2]
        lower = find_zero(λ -> maxloglike1(λ) - τ,(lb[idx_λ],θopt[idx_λ]))
        upper = find_zero(λ -> maxloglike1(λ) - τ,(θopt[idx_λ],ub[idx_λ]), FalsePosition())
        ϵ     = (upper - lower) * 0.01
        P₃ = collect(range(lower+ϵ,upper-ϵ,length=slice))
    end

    # Find the region at each slice
    CR = [[P₃[i], Array{Any,1}(undef,1)] for i = 1:length(P₃)]
    for (i,λ) ∈ enumerate(P₃)

        loglike2 = ξ -> loglike([ξ;λ][invperm([idx_ξ;idx_λ])])
        ξ̂        = optimise(loglike2,θopt[idx_ξ],lb[idx_ξ],ub[idx_ξ])[1]

        # Gradient
        ∇l2      = ξ -> ∇loglike([ξ;λ][invperm([idx_ξ;idx_λ])])[idx_ξ]

        CR[i][2] = confidenceregion2D(ξ̂,loglike2,∇l2;α=α,lb=lb[idx_ξ],ub=ub[idx_ξ],method=method,dt=dt,tmax=tmax,τ=τ)

    end

    CR

end
function confidenceregion3D(θopt,loglike_dv,lb,ub;kwargs...)
    confidenceregion3D(θopt,θ -> loglike_dv(θ)[1], θ -> loglike_dv(θ)[2],lb,ub;kwargs...)
end


function plot_cr3D(CR;slicevar=3,kwargs...)

    idx_λ = slicevar
    idx_ξ = setdiff(1:3,slicevar)

    λ = fill(CR[1][1],size(CR[1][2],1))
    plt = plot([CR[1][2][:,1],CR[1][2][:,2],λ][invperm([idx_ξ;idx_λ])]...,label="";kwargs...)
    plot_cr3D!(plt,CR[2:end];slicevar=slicevar,kwargs...)

end
function plot_cr3D!(plt,CR;slicevar=3,kwargs...)

    idx_λ = slicevar
    idx_ξ = setdiff(1:3,slicevar)

    for i = 1:length(CR)
        λ = fill(CR[i][1],size(CR[i][2],1))
        plot!(plt,[CR[i][2][:,1],CR[i][2][:,2],λ][invperm([idx_ξ;idx_λ])]...,label="";kwargs...)
    end
    return plt

end

function plot_cr3D_proj!(plt,CR;slicevar=3,pα=1.0,kwargs...)

    planes = [xlims(plt)[1],ylims(plt)[2],zlims(plt)[1]]

    idx_λ = slicevar
    idx_ξ = setdiff(1:3,slicevar)

    # Construct matrix of points to plot
    X = vcat([[fill(CR[i][1],size(CR[i][2],1)) CR[i][2]] for i = 1:length(CR)]...)[:,invperm([idx_λ;idx_ξ])]

    # Plot on each axis
    plot!(plt,fill(planes[1],size(X,1)),X[:,2],X[:,3];kwargs...)
    plot!(plt,X[:,1],fill(planes[2],size(X,1)),X[:,3];kwargs...)
    plot!(plt,X[:,1],X[:,2],fill(planes[3],size(X,1));kwargs...)

    # Plot convex hull
    function pts2hull(X)
        pts = [[X[i,1],X[i,2]] for i = 1:size(X,1)]
        hull = convex_hull(pts)
        Y = collect(hcat(hull...)')
        return [Y; Y[[1],:]]
    end
    Y = pts2hull(X[:,[2,3]])
    plot!(plt,fill(planes[1],size(Y,1)),Y[:,1],Y[:,2];kwargs...,α=pα)
    Y = pts2hull(X[:,[1,3]])
    plot!(plt,Y[:,1],fill(planes[2],size(Y,1)),Y[:,2];kwargs...,α=pα)
    Y = pts2hull(X[:,[1,2]])
    plot!(plt,Y[:,1],Y[:,2],fill(planes[3],size(Y,1));kwargs...,α=pα)

end