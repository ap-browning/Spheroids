#=

    total_least_squares.jl

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=#

"""
    lm_orthogonal(X,Y,Z)

Calibrate the "total least squares" line of best fit to 3D data X, Y, Z.
Returns βopt and Lopt, where the line is of the form:
    a,b,ϕ,θ = βopt
    x = [a,b,0] + [cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ)] * t
and where Lopt is the maximised log-likelihood estimate.

Here, the normally distributed error term is taken perpendicular to the line.
"""
function lm_orthogonal(X,Y,Z,s=NaN)

    function orthogonal_error(β,X,Y,Z)
        a,b,ϕ,θ = β

        # Calculate two points on the line
        x₀ = [a,b,0]
        x₁ = x₀ + [cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ)]

        # Calculate errors
        E  = zeros(length(X))
        for i = 1:length(X)
            xᵢ = [X[i],Y[i],Z[i]]
            E[i] = norm((xᵢ - x₀) × (xᵢ - x₁)) / norm(x₁ - x₀)
        end
        E
    end

    # If distribution not provided, fit one first then maximise likelihood
    if (!).(typeof(s) <: Distribution)
        βopt,_ = optimise(β -> -sum(orthogonal_error(β,X,Y,Z)),zeros(4),[-1e6,-1e6,0.0,0.0],[1e6,1e6,π,2π])
        Eopt = orthogonal_error(βopt,X,Y,Z)
        dist = fit_mle(Gamma,Eopt)
        β0   = βopt
    else
        dist = s
        β0   = zeros(4)
    end

    βmle,Lmle = optimise(β -> loglikelihood(dist,orthogonal_error(β,X,Y,Z)),β0,[-1e6,-1e6,0.0,0.0],[1e6,1e6,π,2π])

    return βmle,Lmle,dist

end

"""
    lm_orthogonal_linefcn(β)

Returns x(t) where:
    a,b,ϕ,θ = β
    x(t) = [a,b,0] + [cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ)] * t.

Here, the normally distributed error term is taken perpendicular to the line.
"""
function lm_orthogonal_linefcn(β)
    a,b,ϕ,θ = β
    x₁ = [a,b,0]
    d  = [cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ)]
    t -> x₁ + d * t
end


"""
    lm_orthogonal_test(X,Y,Z,C)

Constructs a likelihood-ratio test to test whether the orthogonal line-of-best-fit
varies with different values of the discrete variable C.

Returns p,β₀,βᵢ where p is the p-value, β₀ is the combined MLE and βᵢ are the MLEs for each value of `sort(unique(C))`.

"""
function lm_orthogonal_test(X,Y,Z,C)

    # Groups
    c = sort(unique(C))

    # Calculate combined MLE (different σ per group) 
    β₀,L₀,d = lm_orthogonal(X,Y,Z)

    # Calculate individual MLEs (store standard deviations)
    βᵢ = Array{Array{Float64,1},1}(undef,length(c))
    Lᵢ = Array{Float64,1}(undef,length(c))
    for i = 1:length(c)
        βᵢ[i],Lᵢ[i],_ = lm_orthogonal(X[C .== c[i]], Y[C .== c[i]], Z[C .== c[i]],d)
    end

    # Seperate values likelihood
    L₁ = sum(Lᵢ)

    # Degrees of freedom in each
    ν₀ = 4
    ν₁ = 4 * length(c)

    # Perform test
    p  = 1.0 - cdf(Chisq(ν₁ - ν₀), 2 * (L₁ - L₀))

    return p,β₀,βᵢ

end