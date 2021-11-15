#=
    Figure 1    
=#

push!(LOAD_PATH,"Module/Greenspan")
push!(LOAD_PATH,"Module/Inference")


using CSV
using DataFrames
using DataFramesMeta
using Plots
using Statistics
using StatsPlots

using Greenspan

include("../FigureDefaults.jl")

#########################################################################
## Fig 1ab: Incucyte (timelapse) data

# Load incucyte data
incucyte = CSV.read("Data/IncucyteData.csv",DataFrame)

# Fig 1a-c 
CellLines = ["983b","793b"]
Densities = [2500,5000,10000]

plts_transient = [plot(title=CellLine) for CellLine in CellLines]
for (i,CellLine) in enumerate(CellLines), (j,Density) in enumerate(Densities)

    data = @subset(incucyte, :CellLine .== CellLine, :InitialCondition .== Density, :Day .≥ [2.0,4.0][i])

    days = unique(data.Day)
    μ    = [mean(filter(isfinite,@subset(data,:Day .== day)[:,:R])) for day in days]
    σ    = [ std(filter(isfinite,@subset(data,:Day .== day)[:,:R])) for day in days]

    plot!(plts_transient[i],days,μ,ribbon=1.96σ,lw=3.0,c=colors[CellLine][Density],xlim=[[2.0,21.0],[4.0,24.0]][i],fα=0.3)

end
fig1ab = plot(plts_transient...,
        link   = :y,
        legend = :none,
        xticks = 0.0:2.0:24.0,
        ylabel = "Radius [μm]",
        xlabel = "Time [d]",
        size   = (800,200)
    )

#########################################################################
## Fig 1cdef

Days = [10,14,18,21]
plts = [plot(title="$t d") for t in Days]
for (i,t) in enumerate(Days), (j,Density) in enumerate(Densities)
    @df @subset(incucyte, :CellLine .== "983b", :InitialCondition .== Density, :Day .== t) density!(plts[i],:R,
        c=colors["983b"][Density],
        frange=0.0,fα=0.3)
end
fig1cdef = plot(plts...,layout=grid(1,4),legend=:none)

#########################################################################
## Fig 1g

    # Solve ODE
    Θ = [0.7,200.0,3.0,0.4,80.0]
    sol = solve_transient_model_all_vars_dim(Θ,20.0)

    T   = range(0.0,20.0,length=200)
    X   = hcat(sol.(T)...)

    R,I,N = X[1,:], X[2,:], X[3,:]

    # Transient solution plot
    fig1g = plot()
    plot!(fig1g, T, R, lw=0, fill=(0,"#99D06C"), label="Cycling")
    plot!(fig1g, T, I, lw=0, fill=(0,"#B077C7"), label="Arrested")
    plot!(fig1g, T, N, lw=0, fill=(0,"#5A6BE5"), label="Necrotic")
    plot!(fig1g, T, R, lw=1.0, c=:black, label="Outer Radius")

    # Phase borders
    t₁ = T[findfirst(R .> Θ[2] * Θ[1])]
    t₂ = T[findfirst(R .> Θ[2])]
    vline!(fig1g,[t₁],c=:black,lw=1.0,ls=:dash,label="")
    vline!(fig1g,[t₂],c=:black,lw=1.0,ls=:dash,label="")

    hoff,voff = 1.0,10.0
    annotate!(fig1g,     0.1 + hoff,10.0 + voff,text("Phase 1",10))
    annotate!(fig1g,t₁ + 0.1 + hoff,10.0 + voff,text("Phase 2",10))
    annotate!(fig1g,t₂ + 0.1 + hoff,10.0 + voff,text("Phase 3",10))

    # Steady state solution plot
    Rinf,_,_ = steady_state(Θ)
    plot!(fig1g,[0.0,20.0],[Rinf,Rinf],ls=:dot,c=:black,lw=2.0,label="Steady State Size")

    plot!(fig1g,
        size   = (700,200),
        legend = :outerright,
        xlim   = (0.0,20.0),
        ylim   = (0.0,320.0),
        xticks = 0:2:20,
        xlabel = "Time [d]",
        ylabel = "Radius [μm]"
    )
    
#########################################################################
## Fig 1
fig1 = plot(fig1ab,fig1cdef,fig1g,layout=grid(3,1),size=(800,600),legend=:none)
#savefig(fig1,"$(@__DIR__)/fig1.svg")
