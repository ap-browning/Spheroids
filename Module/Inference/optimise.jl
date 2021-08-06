#=

    optimise.jl

    Author:     Alexander P. Browning
                ======================
                School of Mathematical Sciences
                Queensland University of Technology
                ======================
                ap.browning@icloud.com
                alexbrowning.me
=#

"""
    optimise(fun,θ₀,lb,ub;kwargs...)

Maximise the function `fun` (vector-valued) within the region bounded by lb × ub

Example:

    fun   = x -> sin(x[1])
    lb,ub = [0.0],[2.0]
    θ₀    = [2 * rand()]
        
    xopt  = optimise(fun,θ₀,lb,ub)

"""
function optimise(fun,θ₀,lb,ub;
        dv = false,
        method = dv ? :LD_LBFGS : :LN_BOBYQA,
    )

    if dv || String(method)[2] == 'D'
        tomax = fun
    else
        tomax = (θ,∂θ) -> fun(θ)
    end

	opt = Opt(method,length(θ₀))
	opt.max_objective = tomax
	opt.lower_bounds = lb       # Lower bound
    opt.upper_bounds = ub       # Upper bound
	opt.local_optimizer = Opt(:LN_BOBYQA, length(θ₀))
	res = optimize(opt,θ₀)
	return res[[2,1]]

end