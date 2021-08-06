#=

    plotting.jl

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=#

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