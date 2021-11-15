#=
    Figure A2

    Total least squares example.
    
=#

using CSV
using DataFrames
using DataFramesMeta
using Plots
using StatsPlots
using Random

using Inference

include("../FigureDefaults.jl")

#########################################################################
## 2D example (y = mx + c)

rng = MersenneTwister(3)

f   = x -> 0.5*x + 2

# Generate data
x   = [1.0,1.9,2.0,3.5,4.6]
y   = f.(x)

n   = length(x)

# Add noise (N(0,0.5))
x += rand(rng,Normal(0.0,1.0),n)
y += rand(rng,Normal(0.0,1.0),n)

# Fit line to data

    # Error function
    function orthogonal_error(β,X,Y)

        m,c = β

        # Calculate errors
        E  = zeros(length(X))
        for i = 1:length(X)
            x₀ = [X[i],Y[i]]
            E[i] = abs(m * X[i] - Y[i] + c) / sqrt(m^2 + 1)
        end
        E

    end

    (m̂,ĉ),_ = optimise(β -> -norm(orthogonal_error(β,x,y)), [2.0,1.0], [0.0,0.0],[10.0,10.0])


# Plot data and line
figA2a = scatter(x,y,c=:black)
        plot!(figA2a,[-2.0,6.0], m̂ * [-2.0,6.0] .+ ĉ, c="#ff375f")

# Plot errors
for i = 1:n

    x₂    = (-(-m̂ * y[i] - x[i]) - m̂ * ĉ) / (m̂^2 + 1)
    y₂    = (m̂ * (x[i] + m̂ * y[i]) + ĉ) / (m̂^2 + 1)

    plot!(figA2a,[x[i],x[i]],[y[i],m̂ * x[i] + ĉ],ls=:dot,c="#007aff")
    plot!(figA2a,[x[i],x₂],[y[i],y₂],c="#007aff")

end

plot!(figA2a,
        legend=:none,
        aspect_ratio=:equal,
        xlabel="x",
        ylabel="y",
        box=:on,
        axis=:all,
        grid=:all,
        xlim=(-1.5,5.5)
    )

#########################################################################
## Error distribution

data = @subset(CSV.read("Data/ConfocalData.csv",DataFrame),:CellLine .== "983b",
            :η .> 0.0, :ϕ .> 0.7, :η .< 0.005 * :R .- 0.85)

# Error function
function orthogonal_error(β,X,Y,Z)
    a,b,ϕ,θ = β

    # Calculate two points on the line
    x₁ = [a,b,0]
    x₂ = x₁ + [cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ)]

    # Calculate errors
    E  = zeros(length(X))
    for i = 1:length(X)
        x₀ = [X[i],Y[i],Z[i]]
        E[i] = norm((x₀ - x₁) × (x₀ - x₂)) / norm(x₂ - x₁)
    end
    E
end
            
# "Orthogonal" linear model
β,L,dist = @df data lm_orthogonal(:R,:ϕ,:η)

# Get errors
E = @df data orthogonal_error(β,:R,:ϕ,:η)

# Plot
figA2b = histogram(E,normalize=:pdf,lw=0.0,c = :Accent_3,label="Error")
plot!(dist,c=:black,label="Fitted Gamma Distribution")
plot!(xlabel="εᵢ",ylabel="Density")

#########################################################################
## Create figure

figA2 = plot(figA2a,figA2b,size=(800,300))
#savefig("$(@__DIR__)/figA2.svg")