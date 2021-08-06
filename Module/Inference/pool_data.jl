#=

    pool_data.jl

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me

=#

"""
    pooled_data(data::DataFrame, groups, vars)

Pool the data in variables `vars` (i.e., subtract the group mean) by `groups`.

"""
function pooled_data(data::DataFrame, groups::Array{Symbol,1}, vars::Array{Symbol,1})

    group_vals  = [unique(data[:,group]) for group ∈ groups]
    filter_vals = collect(Base.product(group_vals...))[:]
    
    pooled_data = Array{Array{Float64,2},1}(undef,length(filter_vals))
    
    # Loop through groups
    for i = 1:length(filter_vals)
    
        ingroup  = .*([data[:,groups[j]] .== filter_vals[i][j] for j = 1:length(groups)]...)
        data_i   = Array(data[ingroup,vars])
        pooled_data[i] = data_i .- mean(data_i,dims=1)
    end
    
    vcat(pooled_data...)

end
function pooled_data(data::DataFrame, groups, vars)
    res = pooled_data(data,
        typeof(groups) <: Array ? groups : [groups],
        typeof(vars)   <: Array ? vars   : [vars]
    )
    if typeof(vars) == Symbol
        return res[1]
    else
        return res
    end
end

pooled_std(data::DataFrame, groups, vars) = std(pooled_data(data,groups,vars),dims=1)[:]
pooled_cov(data::DataFrame, groups, vars) = cov(pooled_data(data,groups,vars))
pooled_cor(data::DataFrame, groups, vars) = cor(pooled_data(data,groups,vars))