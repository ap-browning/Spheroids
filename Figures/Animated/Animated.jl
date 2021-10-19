#=
    Animation of Greenspan model through time.    
=#

using CSV
using DataFrames
using DataFramesMeta
using Plots
using Statistics
using StatsPlots

using Greenspan

# Parameters in our "reduced" parameterisation
Q = 0.7
R₁ = 200.0
γ = 3.0
s = 0.4
R₀ = 80.0

# Parameters
λ = γ * s * 3
ω₂ = 100.0
ω₁ = 20.0
k = 1000.0
A = 12.0
β₁ = 300.0
κ = 500.0 
P = A / (k * (ω₂ - ω₁)) * β₁ * κ / Q^2

# Solve model
sol = solve_transient_model_all_vars_dim([Q,R₁,γ,s,R₀],20.0)

# Function to plot circle
function circle(r)
    θ = range(0,2π,length=100)
    r*sin.(θ),r*cos.(θ)
end

# Setup grid
x = range(-350.0,350.0,length=100);
y = range(-350.0,350.0,length=100);
R = [sqrt(xx^2 + yy^2) for xx in x, yy in y];

# Oxygen and Nutrients as a function of r
function ωfun(r,Ro,Ri,Rn)
    if r > Ro
        return ω₂
    elseif Rn ≤ r ≤ Ro
        return ω₂ - A / 6k * (Ro^2 - r^2) + A * Rn^3 / 3k * (1 / r - 1 / Ro)
    else
        return ω₁
    end
end
function βfun(r,Ro,Ri,Rn)
    if Rn ≤ r ≤ Ro
        return P / 3κ * (0.5 * (Ro^2 - r^2) - Rn^3 * (1 / r - 1 / Ro))
    elseif r < Rn
        return βfun(Rn,Ro,Ri,Rn)
    else
        return 0.0
    end
end

# Produce animation
anim = @animate for tp in range(0.0,20.0,step=0.1)
	
    Ro,Ri,Rn = sol(tp)
    p1 = plot(circle(Ro),c="#A6D272",frange=0.0,label="Cycling")
    plot!(p1,circle(Ri),c="#AC79B5",frange=0.0,label="Inhibited")
    plot!(p1,circle(Rn),c="#5C6AB2",frange=0.0,label="Necrotic")
    plot!(p1,aspect_ratio=:equal,xlim=(-350.0,350.0),ylim=(-350.0,350.0))
    plot!(p1,title="$tp d",box=:on)

    r = range(0.0,350.0,length=100)
        
    # Nutrients
    ω = ωfun.(R,Ro,Ri,Rn)
    p2 = heatmap(x,y,ω',clim=(0.0,100.0),aspect_ratio=:equal,xlim=(-350.0,350.0),ylim=(-350.0,350.0),title="Nutrient",colorbar=false)
    plot!(p2,box=:on,xticks=[],yticks=[],legend=:none)

    # Waste
    β = βfun.(R,Ro,Ri,Rn)
    p3 = heatmap(x,y,β',clim=(0.0,1000.0),aspect_ratio=:equal,xlim=(-350.0,350.0),ylim=(-350.0,350.0),title="Waste",c=:acton,colorbar=false)
    plot!(p3,box=:on,xticks=[],yticks=[],legend=:none)

    plot!(p2,circle(Ro),c=:black)
    plot!(p3,circle(Ro),c=:black)
    if Ri > 0.0
        plot!(p2,circle(Ri),c=:black)
        plot!(p3,circle(Ri),c=:black)
    end
    if Rn > 0.0
        plot!(p2,circle(Rn),c=:black)
        plot!(p3,circle(Rn),c=:black)
    end

    p4 = plot([0.0,Rn],[110.0,110.0],frange=-10.0,lw=0.0,fc="#5C6AB2",fα=0.5)
    plot!(p4,[Rn,Ri],[110.0,110.0],frange=-10.0,lw=0.0,fc="#AC79B5",fα=0.5)
    plot!(p4,[Ri,Ro],[110.0,110.0],frange=-10.0,lw=0.0,fc="#A6D272",fα=0.5)
    hline!(p4,[ω₁],ls=:dash,lw=2.0,c=:grey)
    plot!(p4,r,ωfun.(r),ylim=(-1.0,101.0),c=:black,label="",lw=2.0,xlim=(0.0,350.0),legend=:none,xlabel="Radius",ylabel="Nutrient")
    annotate!(p4,300.0,25.0,"ωcrit",fontfamily="Helvetica")

    p5 = plot([0.0,Rn],[1100.0,1100.0],frange=-10.0,lw=0.0,fc="#5C6AB2",fα=0.5)
    plot!(p5,[Rn,Ri],[1100.0,1100.0],frange=-10.0,lw=0.0,fc="#AC79B5",fα=0.5)
    plot!(p5,[Ri,Ro],[1100.0,1100.0],frange=-10.0,lw=0.0,fc="#A6D272",fα=0.5)
    hline!(p5,[β₁],ls=:dash,lw=2.0,c=:grey)
    plot!(p5,r,βfun.(r),ylim=(-10.0,1010.0),c=:black,label="",lw=2.0,xlim=(0.0,350.0),legend=:none,xlabel="Radius",ylabel="Waste")
    annotate!(p5,300.0,360.0,"βcrit",fontfamily="Helvetica")

    plot(p1,p2,p3,p4,p5,layout=@layout([grid(1,3);grid(1,2)]),size=(800.0,600.0),fontfamily="Helvetica")

end

mp4(anim,"Figures/Animated/Animation.mp4",fps=10)