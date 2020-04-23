"""
	compute_fitness!(traj::Array{<:FrequencyTraj,1}, fp::FluPop, ftype, trajfield; strainfield=:lbi, shift=:mean)
"""
function compute_fitness!(traj::Array{<:FrequencyTraj,1}, fp::FluPop, ftype; trajfield=Symbol(ftype,:_fitness), strainfield=:lbi, shift=:mean)
	if trajfield == ftype
		@error "Cannot replace $trajfield field when computing fitness."
	end
	_tf = (ftype == :strains && trajfield == Symbol(ftype,:_fitness)) ? Symbol(strainfield,:_fitness) : trajfield
	# 
	for t in traj
		t.data[_tf] = zeros(Float64, length(t))
	end
	# 
	if ftype == :strains
		compute_strains_fitness!(traj, fp, strainfield, trajfield=_tf, shift=shift)
	elseif ftype == :date
		for t in traj
			compute_date_fitness!(t, _tf)
		end
	elseif ftype == :region
		for t in traj
			ref_regions = Dict(datebin_to_date(db)=>get_regions(s) for (db,s) in fp.datebin)
			get_regions!(t, fp)
			compute_region_fitness!(t, ref_regions, _tf)
		end
	elseif ftype == :treespread
		for t in traj
			compute_treespread_fitness!(t, _tf, score = :squaredfreqs)
		end
	end
end

"""
	compute_date_fitness!(traj::FrequencyTraj)

```
traj.fitness[i] = traj.t[i] - traj.t[1]
```
"""
function compute_date_fitness!(traj::FrequencyTraj, field)
	for (i,t) in enumerate(traj.t)
		traj.data[field][i] = t.value - traj.t[1].value
	end
end


"""
	compute_strains_fitness!(traj::FrequencyTraj, fp::FluPop, strainfield=:lbi, trajfield=strainfield; shift = :mean)
	compute_strains_fitness!(traj::Array{<:FrequencyTraj,1}, fp::FluPop, strainfield=:lbi, trajfield=strainfield; shift = :mean)
"""
function compute_strains_fitness!(traj::FrequencyTraj, fp::FluPop, strainfield=:lbi; trajfield=Symbol(strainfield,:_fitness), shift = :mean)
	for (i, strains) in enumerate(traj.strains)
		current_date = traj.date + traj.t[i]
		db = find_datebin(current_date, fp)
		# Strains in trajectory
		Ntraj = count(x->!ismissing(fp.strains[x][strainfield]), strains)
		Ftraj = mapreduce(x->ismissing(fp.strains[x][strainfield]) ? 0. : fp.strains[x][strainfield], +, strains, init=0.) 
		# All strains in the datebin
		N = count(x->!ismissing(x[strainfield]), fp.datebin[db])
		F = mapreduce(x->ismissing(x[strainfield]) ? 0. : x[strainfield], +, fp.datebin[db], init=0.) 
		if ismissing(F) || ismissing(Ftraj)
			traj.data[trajfield][i] = missing # This would be the proper way to do it. 
		elseif Ntraj == 0 || Ntraj == N
			traj.data[trajfield][i] = 0.
		elseif N > Ntraj
			if shift==:exclusive_mean
				traj.data[trajfield][i] = Ftraj/Ntraj  -  (F - Ftraj)/(N - Ntraj)
			elseif shift==:mean
				traj.data[trajfield][i] = Ftraj/Ntraj - F/N
			elseif shift==:none
				traj.data[trajfield][i] = Ftraj/Ntraj
			else
				@error "Possible values of `shift` are `:none`, `:mean`, `:exclusive_mean`."
			end
		else
			@error "Ntraj $(Ntraj) > N $N - More strains strains in trajectory than in datebin."
		end
	end
end
function compute_strains_fitness!(traj::Array{<:FrequencyTraj,1}, fp::FluPop, strainfield=:lbi; trajfield=Symbol(strainfield,:_fitness), shift = :mean)
	for t in traj
		compute_strains_fitness!(t, fp, strainfield, trajfield=trajfield, shift=shift)
	end
end


"""
	 compute_treespread_fitness!(traj::FrequencyTraj; score = :squaredfreqs)

Default : if `traj.treespread` is empty, the trajectory is considered as completely localized in the tree for the corresponding timebin (*i.e.* of lowest fitness).
"""
function compute_treespread_fitness!(traj::FrequencyTraj, field; score = :squaredfreqs)
	for (i,ts) in enumerate(traj.data[:treespread])
		x = ts / sum(ts)
		if score == :squaredfreqs
			traj.data[field][i] = isempty(x) ? 0. : 1 - sum(x.^2)
		elseif score == :entropy
			traj.data[field][i] = isempty(x) ? 0. : StatsBase.entropy(x)
		end
	end
end

"""
	compute_region_fitness!(traj::FrequencyTraj, field)
"""
function compute_region_fitness!(traj::FrequencyTraj, ref_regions::Dict, field)
	for (i,r) in enumerate(traj.data[:regions])
		traj.data[field][i] = (isempty(ref_regions) || isempty(r)) ? 0. : (compute_region_entropy(r))# - compute_region_entropy(ref_regions[traj.date + traj.t[i]]))
	end
end
function compute_region_entropy(regionspread; scoretype=:entropy)
	x = collect(values(regionspread))
	!isempty(x) && (x /= sum(x))
	if scoretype == :squaredfreqs
		return isempty(x) ? 0. : 1 - sum(x.^2)
	elseif scoretype == :entropy
		return StatsBase.entropy(x)
	else
		@error "`scoretype` can be `:squaredfreqs` or `:entropy`"
	end
end