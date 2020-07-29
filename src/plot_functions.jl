using Distributions, Plots
using FluPredictibility.BioTools

export plot_all_trajectories
export pfix, pfix_v_freq, pfix_v_freq_positivederivative, fitness_plot, meanfreq


function plot_all_trajectories(ph; yearticks=false, label=true, lw=2.5)
	p = plot(size=(900,600))
	for z in ph
		X,Y,tmp,α = frequency_series(z)
		if yearticks
			X = [year(x) + month(x) /12. for x in X]
		end
		for a in 1:size(Y,2)
			plot!(p, X, Y[:,a], label=label ? "$(α[a])" : "", linewidth=lw)
		end
	end
	plot!(p, ylim=(-0.05,1.05))
	return p
end

function plot_all_trajectories(trajectories::Array{FrequencyTraj{A},1}; label=false, lw = 1) where A
	p = plot(size=(900,600))
	for traj in trajectories
		X = traj.date .+ traj.t 
		Y = traj.freq
		lab = "$(traj.i) - $(traj.val)"
		if traj.fixation == :fixed
			st = (lw, :blue, 0.5)
		elseif traj.fixation == :lost
			st = (lw, :red, 0.5)
		else
			st = (lw, :black, :dashdot, 0.5)
		end
		plot!(p, X, Y, label = label ? lab : "", line=st)
	end
	plot!(p, frame=:box)
	return p
end


function trajectory_freqbin(traj, alphabins)
	freqtraj_cf = Dict()
	for (α,dα) in alphabins
	    freqtraj_cf[α] = frequency_condition(traj, α, dα=dα)
	end
	return freqtraj_cf
end

function bernoulli_estimator(x,n)
	naivemean = x/n
	estmean = (x+1/2) / (n+1)
	α = x + 1/2
	β = n - x + 1/2
	P = Distributions.Beta(α,β)
	# 
	lowerbound = 0.
	while cdf(P, lowerbound) < 0.05
		lowerbound += 0.01
	end
	# 
	higherbound = 1.
	while cdf(P, higherbound) > 0.95
		higherbound -= 0.01
	end
	# 
	return estmean, estmean-lowerbound, higherbound-estmean
end

# pfix(traj) = sum([x.fixation==:fixed for x in traj]) / (sum([x.fixation==:fixed for x in traj]) + sum([x.fixation==:lost for x in traj]));
meanfreq(traj) = mean(t.freq[t.index[:active]] for t in traj)
"""
"""
function pfix(trajectories)
	n = length(trajectories)
	x = count(t->t.fixation==:fixed, trajectories)
	yf, errdown, errup = bernoulli_estimator.(x,n)
	# 
	xf = meanfreq(trajectories)
	err = (errdown, errup)
	return xf, yf, err
end

"""
"""
function pfix_v_freq(ph, alphabins; v=false)
	trajectories = all_trajectories(ph, keep_unfinished=false)
	trajectories = previous_state_condition(trajectories, :lost)
	# Binning by frequency
	traj_fb = sort(OrderedDict(trajectory_freqbin(trajectories, alphabins)));
	# Keeping only trajectories that have a frequency backed by 20 strains at the time where it is binned. 
	for (k,v) in traj_fb
	    traj_fb[k] = population_size_condition(v, 20, mode=:active)
	end
	#
	v && println("Based on a total of $(sum(length(v) for v in values(traj_fb))) trajectories")
	# 
	n = [length(traj_fb[x]) for x in keys(traj_fb)] # For error bars
	x = [count(t->t.fixation==:fixed, traj_fb[x]) for x in keys(traj_fb)]
	out = bernoulli_estimator.(x,n)
	# 
	xf = [meanfreq(traj_fb[x]) for x in keys(traj_fb)]
	yf = [x[1] for x in out]; 
	# yf = [count(t->t.fixation==:fixed, traj_fb[x])/length(traj_fb[x]) for x in keys(traj_fb)]
	errup = [x[3] for x in out]; errdown = [x[2] for x in out]
	# 
	err = (errdown, errup)
	# return sort(collect(keys(traj_fb))), yf, err
	# println(xf)
	return xf, yf, err
end

"""
	fitness_plot(trajectories, field, alphabins; verbose=false, dq = 0., alldat=false)

Trajectories higher or lower than median fitness.
"""
function fitness_plot(trajectories, field, alphabins; verbose=false, dq = 0., alldat=false)
    traj_fb = trajectory_freqbin(trajectories, alphabins);
    # 
    for (k,v) in traj_fb
        traj_fb[k] = population_size_condition(v, 20, mode=:active)
    end   
    # High and low fitnesses
    high_fit = Dict()
    low_fit = Dict()
    medfit = Dict()
    for (k,v) in traj_fb
        fvalues = [x.data[field][x.index[:active]] for x in v]
        medfit[k] = median(fvalues)
        # println(fvalues)
        verbose && println("Frequency $k -- median fitness $(medfit[k])")
        high_fit[k] = v[findall(x->x.data[field][x.index[:active]] > quantile(fvalues, 0.5 + dq), v)]
        low_fit[k] = v[findall(x->x.data[field][x.index[:active]] <= quantile(fvalues, 0.5 - dq), v)]
    end
    # Arrays
    dat = vcat([reshape(collect(pfix(traj_fb[x[1]])), 1, 3) for x in alphabins]...)
    dat_low = vcat([reshape(collect(pfix(low_fit[x[1]])), 1, 3) for x in alphabins]...)
    dat_high = vcat([reshape(collect(pfix(high_fit[x[1]])), 1, 3) for x in alphabins]...)
	#
    alldat && return dat, dat_low, dat_high, traj_fb, high_fit, low_fit, medfit
    return dat, dat_low, dat_high
end
"""
	fitness_plot(trajectories, field, alphabins, fitval; verbose=false alldat=false)

Trajectories higher or lower than `fitval` for field `field`. 
"""
function fitness_plot(trajectories, field, alphabins, fitval; verbose=false, alldat=false)
    traj_fb = trajectory_freqbin(trajectories, alphabins);
    # 
    for (k,v) in traj_fb
        traj_fb[k] = population_size_condition(v, 20, mode=:active)
    end   
    # High and low fitnesses
    high_fit = Dict()
    low_fit = Dict()
    medfit = Dict()
    for (k,v) in traj_fb
        high_fit[k] = v[findall(x->x.data[field][x.index[:active]] > fitval, v)]
        low_fit[k] = v[findall(x->x.data[field][x.index[:active]] <= fitval, v)]
    end
    # Arrays
    dat = vcat([reshape(collect(pfix(traj_fb[x[1]])), 1, 3) for x in alphabins]...)
    dat_low = vcat([reshape(collect(pfix(low_fit[x[1]])), 1, 3) for x in alphabins]...)
    dat_high = vcat([reshape(collect(pfix(high_fit[x[1]])), 1, 3) for x in alphabins]...)
	#
    alldat && return dat, dat_low, dat_high, traj_fb, high_fit, low_fit, medfit
    return dat, dat_low, dat_high
end

"""
"""
function pfix_v_freq_positivederivative(ph, alphabins)
	trajectories = all_trajectories(ph, keep_unfinished=false)
	trajectories = previous_state_condition(trajectories, :lost)
	# Binning by frequency
	traj_fb = sort(OrderedDict(trajectory_freqbin(trajectories, alphabins)));
	# Keeping only trajectories that have a frequency backed by 50 strains at the time where it is binned. 
	for (k,v) in traj_fb
	    traj_fb[k] = population_size_condition(v, 20, mode=:active)
	    traj_fb[k] = derivative_condition(traj_fb[k])
	end
	# 
	n = [length(traj_fb[x]) for x in keys(traj_fb)] # For error bars
	x = [count(t->t.fixation==:fixed, traj_fb[x]) for x in keys(traj_fb)]
	out = bernoulli_estimator.(x,n)
	# 
	xf = [meanfreq(traj_fb[x]) for x in keys(traj_fb)]
	yf = [x[1] for x in out]; 
	errup = [x[3] for x in out]; errdown = [x[2] for x in out]
	# 
	err = (errdown, errup)
	return xf, yf, err
end

function mean_trajectory()

end	
