#=

    steady_state.jl

    Contains functions required to solve the steady state Greenspan model.

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


    Parameterisation:
        R ∈ (0,∞) The outer radius
        ϕ ∈ (0,1) The ratio of inhibited to outer radius
        η ∈ (0,1) The ratio of necrotic to outer radius

        x = [R,ϕ,η] at any time
        X = [R,ϕ,η] x specifically at steady state

        ρ ∈ (0,1) An intermediate parameter, equal to η / ϕ

    Parameters:
        Q  ∈ (0,1)
        γ  ∈ (0,∞)
        R₁ ∈ (0,∞)

        θ = [Q,R₁,γ]  A Float64 array of all parameters

    We note that
        [ϕ,η] = f(R;θ)
    at any time.

=#

"""
    _check_θ(θ)

Verifies that parameter inputs are valid.
"""
function _check_θ(Θ)

    Q,R₁,γ  = Θ[1:3]

    (Q  ≥ 0 && Q ≤ 1) || error("Invalid value Q = $Q.")
    (γ  > 0)          || error("Invalid value γ = $γ.")
    (R₁ ≥ 0)          || error("Invalid value R₁ = $R₁.")

    return Θ[1:3]

end

"""
    coefs(θ)

Computes the coefficients of the polynomial for ρ.

Example: 
    f(ρ) = Polynomial(coefs([0.8,0.5,150.0]))(ρ)
"""
function coefs(θ)

    Q,γ = θ[[1,3]]

    c = [
        3Q^2 - 3Q^4 + Q^6,
        0,
        -9Q^2,
        18Q^2 - 18Q^4 + 6Q^6 - 2γ+ 9Q^2*γ - 9Q^4*γ + 3Q^6*γ,
        27Q^4,
        -36Q^2 - 9Q^2*γ,
        36Q^2 - 36Q^4 - 15Q^6 - 6γ + 36Q^2*γ - 36Q^4*γ + 12Q^6*γ - 3γ^2 + 9Q^2*γ^2 - 9Q^4*γ^2 + 3Q^6*γ^2,
        54Q^4 + 27Q^4*γ,
        -36Q^2 - 36Q^2*γ,
        24Q^2 - 24Q^4 + 8Q^6 + 36Q^2*γ - 36Q^4*γ - 15Q^6*γ - 6γ^2 + 18Q^2*γ^2 - 18Q^4*γ^2 + 6Q^6*γ^2 - γ^3 + 3Q^2*γ^3 - 3Q^4*γ^3 + Q^6*γ^3,
        54Q^4*γ,
        -36Q^2*γ,
        8γ
    ]

end

"""
    steady_state(θ)

Computes the steady state, M.

Example: 
    M = steady_state([0.8,0.5,150.0])
"""
function steady_state(Θ)

    θ = _check_θ(Θ)

    # Solve polynomial for ρ
    ρ     = find_zero(Polynomial(coefs(θ)),(0,1))

    # Calculate ϕ
    fϕ    = (ρ,Q,R₁,γ) -> (1.0 + γ*ρ^3)^(-1/3)
    ϕ     = fϕ(ρ,θ...)

    # Calculate R
    fR    = (ρ,ϕ,Q,R₁,γ) ->  R₁ * (1.0 - 3ρ^2*ϕ^2 + 2ρ^3*ϕ^3)^(-1/2)
    R     = fR(ρ,ϕ,θ...)

    # Calculate η
    η     = ϕ * ρ

    return [R,ϕ,η]

end


"""
    steady_state_dv(θ)

Computes the steady state, M, and the Jacobian, J.

Example: 
    M,J = steady_state([0.8,0.5,150.0])
"""
function steady_state_dv(Θ)

    θ = _check_θ(Θ)

    # Polynomial coefficients
    c     = coefs(θ)

    # Solve polynomial for ρ ,Roots.FalsePosition()
    ρ     = find_zero(Polynomial(coefs(θ)),(0,1))

    # Calculate ϕ
    fϕ    = (ρ,Q,R₁,γ) -> (1.0 + γ*ρ^3)^(-1/3)
    ϕ     = fϕ(ρ,θ...)

    # Calculate R
    fR    = (ρ,ϕ,Q,R₁,γ) ->  R₁ * (1.0 - 3ρ^2*ϕ^2 + 2ρ^3*ϕ^3)^(-1/2)
    R     = fR(ρ,ϕ,θ...)

    # Calculate η
    η     = ϕ * ρ

    # Calculate Jacobian ∂(R,ϕ,ρ)/∂(Q,γ,R₁)

        # ∂cᵢ/∂θ
        ∂c    = ForwardDiff.jacobian(coefs,θ)

        # ∂ρ/∂cᵢ
        i     = 0:12
        ∂ρ∂cᵢ = -ρ.^i ./ sum(c .* i .* ρ.^(i.-1))

        # Calculate ∇ρ
        ∇ρ    = (∂ρ∂cᵢ' * ∂c)[:]

        # Calculate ∇ϕ
        ∂ϕ∂ρ  = ForwardDiff.derivative(ρ -> fϕ(ρ,θ...),ρ)
        ∂ϕ∂θ  = ForwardDiff.gradient(θ -> fϕ(ρ,θ...),θ)
        ∇ϕ    = ∂ϕ∂ρ * ∇ρ + ∂ϕ∂θ

        # Calculate ∇R
        fR    = (ρ,ϕ,Q,R₁,γ) ->  R₁ * (1.0 - 3ρ^2*ϕ^2 + 2ρ^3*ϕ^3)^(-1/2)
        ∂R∂ρ  = ForwardDiff.derivative(ρ -> fR(ρ,ϕ,θ...),ρ)
        ∂R∂ϕ  = ForwardDiff.derivative(ϕ -> fR(ρ,ϕ,θ...),ϕ)
        ∂R∂θ  = ForwardDiff.gradient(θ -> fR(ρ,ϕ,θ...),θ)
        ∇R    = ∂R∂ρ * ∇ρ + ∂R∂ϕ * ∇ϕ + ∂R∂θ

        # Calculate ∇η
        ∂η∂ρ  = ϕ
        ∂η∂ϕ  = ρ
        ∇η    = ∂η∂ρ * ∇ρ + ∂η∂ϕ * ∇ϕ

        J     = [∇R ∇ϕ ∇η]'
        
    return [R,ϕ,η],J

end