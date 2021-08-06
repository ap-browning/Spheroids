#=

    transient_model.jl

    Contains functions required to solve the transient Greenspan model.

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

    Parameters:
        Q  ∈ (0,1)
        γ  ∈ (0,∞)
        R₁ ∈ (0,∞)

        θ = [Q,R₁,γ]  A Float64 array of all parameters
        Θ = [θ;s,R₀]  Include s, the transient variable

    We note that
        [ϕ,η] = f(R;θ)
    at any time.

=#

"""
    _inner_vars(R,θ)

Calculate the inner variables (ϕ,η) for a fixed R given parameters θ

Example:
    ϕ,η = _inner_vars(200.0,[0.8,0.5,150.0])
"""
function _inner_vars(R::Number,θ)

    Q,R₁,γ = _check_θ(θ[1:3])

    # Phase 1
    if R ≤ Q * R₁
        ϕ,η = 0.0,0.0

    # Phase 2
    elseif R ≤ R₁
        ϕ   = sqrt(1 - (Q*R₁)^2 / R^2)
        η   = 0.0

    # Phase 3
    elseif R > R₁
        η   = find_zero(Polynomial([R^2 - R₁^2, 0.0, -3R^2, 2R^2]),[0.0,1.0])
        ϕ   = find_zero(Polynomial([2*η^3*R^2,Q^2*R₁^2 - R^2*(1 + 2*η^3),0.0,R^2]),[η,1.0])

    else
        error("Invalid phase.")
    end

    return [ϕ,η]

end
function _inner_vars(R::Vector,θ)

    x = hcat([_inner_vars(r,θ) for r = R]...)
    Φ = x[1,:]
    Η = x[2,:]

    return Φ,Η

end

"""
    _ode(R,Θ,t)

Calculate the right hand side of the Greenspan ODE model.

Example:
    dR = _ode(100.0,[0.8,0.5,150.0,1.0],0.0)
"""
function _ode(R,Θ,t)

    # Parameters
    Q,R₁,γ,s = Θ[1:4] # Note capital Θ

    # Inner variables
    ϕ,η = _inner_vars(R,Θ)

    # Calculate RHS
    dR = s/3 * R * (1 - ϕ^3 - γ*η^3)

end


"""
    solve_transient_model(Θ,tmax)

Solve the Greenspan model from t = 0.0 to t = tmax with parameters
    Θ = [Q,R₁,γ,s,R₀]

Example:
    R = solve_transient_model([0.8,0.5,150.0,1.0],20.0)
    plot(R)
"""
function solve_transient_model(Θ,tmax)
    solve(ODEProblem(_ode,Θ[end],(0.0,tmax),Θ[1:4]))
end

"""
    solve_transient_model_all_vars(Θ,tmax)

Solve the Greenspan model from t = 0.0 to t = tmax with parameters
    Θ = [Q,R₁,γ,s,R₀]
and return a function F, where F(t) = [R(t),ϕ(t),η(t)]

Example:
    F = solve_transient_model_all_vars([0.8,0.5,150.0,1.0],20.0)
"""
function solve_transient_model_all_vars(ode_sol,Θ)

    function all_vars(t)
        R = ode_sol(t)
        ϕ,η = _inner_vars(R,Θ)
        [R,ϕ,η]
    end
    all_vars

end
solve_transient_model_all_vars(Θ::Vector,tmax) = solve_transient_model_all_vars(solve_transient_model(Θ,tmax),Θ)

"""
    solve_transient_model_all_vars_dim(Θ,tmax)

Solve the Greenspan model from t = 0.0 to t = tmax with parameters
    Θ = [Q,R₁,γ,s,R₀]
and return a function F, where F(t) = [R(t),R(t) * ϕ(t),R(t) * η(t)],
i.e., the dimensional parameters.

Example:
    F = solve_transient_model_all_vars_dim([0.8,0.5,150.0,1.0],20.0)
"""
function solve_transient_model_all_vars_dim(ode_sol,Θ)
    function all_vars(t)
        R = ode_sol(t)
        ϕ,η = _inner_vars(R,Θ)
        [R,R*ϕ,R*η]
    end
    all_vars
end
solve_transient_model_all_vars_dim(Θ::Vector,tmax) = solve_transient_model_all_vars_dim(solve_transient_model(Θ,tmax),Θ)