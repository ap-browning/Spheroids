#=

    structure_model.jl

    Contains functions required to solve the Greenspan structure.

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

    Parameters:
        Q  ∈ (0,1)
        γ  ∈ (0,∞)         
        R₁ ∈ (0,∞)

        θ = [Q,R₁,γ]  A Float64 array of all parameters
        ξ = [Q,R₁]    A Float64 array of parameters in f_struct

    We note that
        [ϕ,η] = f(R;θ)
    at any time.

=#

"""
    structure_model(ξ,R)

Solve the structure model for a given ξ = [Q,R₁], given R.

Example: 
    ϕ,η = structure_model([0.8,150.0],200.0)
"""
function structure_model(ξ::Vector,R::Number)
    _inner_vars(R,[ξ;1.0])
end
function structure_model(ξ::Vector,R::Vector)
    call = Rᵢ -> _inner_vars(Rᵢ,[ξ;1.0])
    call.(R)
end

"""
    structure_model_dv(ξ,R)

Solve the structure model for a given ξ = [Q,R₁], given R.
Additionally return the model Jacobian.

Example: 
    M,J = structure_model_dv([0.8,150.0],200.0)
"""
function structure_model_dv(ξ::Vector,R::Number)

    Q,R₁ = ξ

    # Phase 1
    if R ≤ Q * R₁
        return zeros(2), zeros(2,2)

    # Phase 2
    elseif R ≤ R₁
        ϕ  = sqrt(1 - (Q*R₁)^2 / R^2)
        ∇ϕ = -Q * R₁ / (ϕ * R^2) * [R₁,Q]
        J  = [∇ϕ zeros(2)]'
        return [ϕ,0.0],J
        
    # Phase 3
    elseif R > R₁

        function coefs_η(ξ)
            Q,R₁ = ξ
            [R^2 - R₁^2, 0.0, -3R^2, 2R^2]
        end
        function coefs_ϕ(ξ,η)
            Q,R₁ = ξ
            [2*η^3*R^2,Q^2*R₁^2 - R^2*(1 + 2*η^3),0.0,R^2]
        end

        c_η = coefs_η(ξ)
        η   = find_zero(Polynomial(c_η),[0.0,1.0])

        c_ϕ = coefs_ϕ(ξ,η)
        ϕ   = find_zero(Polynomial(c_ϕ),[η,1.0])

        ∂c_η∂ξ = ForwardDiff.jacobian(coefs_η,ξ)
        ∂c_ϕ∂ξ = ForwardDiff.jacobian(ξ -> coefs_ϕ(ξ,η),ξ)
        ∂c_ϕ∂η = ForwardDiff.derivative(η -> coefs_ϕ(ξ,η),η)

        i = 0:3
        ∂η∂c_ηᵢ = -η .^ i / sum(c_η .* i .* η.^(i.-1.0)) 
        ∂ϕ∂c_ϕᵢ = -ϕ .^ i / sum(c_ϕ .* i .* ϕ.^(i.-1.0)) 

        ∇η = (∂η∂c_ηᵢ' * ∂c_η∂ξ)[:]
        ∇ϕ = (∂ϕ∂c_ϕᵢ' * ∂c_ϕ∂ξ)[:] + (∂ϕ∂c_ϕᵢ' * ∂c_ϕ∂η) * ∇η

        J = [∇ϕ ∇η]'

        return [ϕ,η],J

        η   = find_zero(Polynomial([R^2 - R₁^2, 0.0, -3R^2, 2R^2]),[0.0,1.0])
        ϕ   = find_zero(Polynomial([2*η^3*R^2,Q^2*R₁^2 - R^2*(1 + 2*η^3),0.0,R^2]),[η,1.0])

    else
        error("Invalid phase.")
    end

end