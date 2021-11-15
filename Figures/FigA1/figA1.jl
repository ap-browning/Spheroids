#=
    Figure A1

    (a) Number of roots to the polynomial describing steady-state solution
    (b) Semi-analytical solution vs long-term transient solution

=#

using Plots
using Polynomials

using Greenspan

include("../FigureDefaults.jl")

#########################################################################
## (a) Plot number of roots of f(ρ; Q, γ)

Qrange = range(-0.5,1.5,length=200)
γrange = range(-5.0,20.0,length=200)
Nroots = zeros(length(Qrange),length(γrange))

for (i,Q) ∈ enumerate(Qrange), (j,γ) ∈ enumerate(γrange)

    # Obtain all 12 roots
    all_roots = roots(Polynomial(coefs([Q,100.0,γ])))

    # Count number of real roots between 0.0 and 1.0
    Nroots[i,j] = count((imag(all_roots) .== 0) .& (0.0 .≤ real(all_roots) .≤ 1.0))

end

figA1a = plot(Qrange,γrange,Nroots',st=:contourf,lw=0.0,clim=(0.0,3.0),
    ylim=(-5.0,20.0),xlim=(-0.5,1.5),axis=:on,
    xlabel="Q [-]", ylabel="Gamma [-]", colorbar_title="Number of solutions",
    c = :YlGnBu_3
)
plot!(figA1a,[0.0,1.0,1.0,0.0,0.0],[0.0,0.0,30.0,30.0,0.0],c=:black,legend=:none)


#########################################################################
## (b) Plot solution to transient model, and steady state solution

# Parameters (Θ = [Q,γ,R₁,s,R₀])
Θ = [0.8,150.0,1.0,1.0,100.0]

# Obtain transient solution
sol = solve_transient_model_all_vars_dim(Θ,50.0)

# Obtain steady state solution, and scale up by R
SS  = steady_state(Θ); SS_dim = [SS[1]; SS[2:3] .* SS[1]]

# Plot
T   = range(0.0,50.0,length=200)
figA1b = plot(T,hcat(sol.(T)...)',labels=["R(t)" "ϕ(t)" "η(t)"])
    hline!(figA1b,[SS_dim[1]],c=:black, ls=:dash, label="R̄")
    hline!(figA1b,[SS_dim[2]],c=:black, ls=:dashdot, label="R̄ϕ̄")
    hline!(figA1b,[SS_dim[3]],c=:black, ls=:dot, label="R̄η̄")
    plot!(figA1b,
        xlabel = "Time [d]",
        ylabel = "Radius [um]",
        legend = :bottomright
    )

#########################################################################
## Produce figure
figA1 = plot(figA1a,figA1b,size=(800,220),layout=@layout([a{0.4w} b]))
#savefig("$(@__DIR__)/figA1.svg")
