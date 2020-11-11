"""
	all_trajectories(posh::PosEvo; fixed_thr = 0.95, lost_thr = 0.05, keep_unfinished=false)

Find all active trajectories for given position `posh`. 
An active trajectory is defined by being above fixing threshold (resp. below absence threshold) for 2 timebins in a row, and below (resp. above) at the next timebin. 
Those can then be filtered for given characteristics, such as being seen in some frequency bin, being absent or fixed in the past, etc...  

## Technical note
`start_indices`: `Y[i-1,a]>=fixed_thr && Y[i,a]>=fixed_thr` (*i.e.* mutation was fixed) AND `Y[i+1,a]<fixed_thr, *i.e.* it is now not fixed. This and the reverse for `lost_thr`.   
`end_indices`: `Y[i,a] >= fixed_thr && Y[i+1,a] >= fixed_thr` (*i.e.* mutation is now fixed) AND `Y[i-1,a] < fixed_thr` (*i.e.* it was not fixed before). 
"""
function all_trajectories(posh::PosEvo{A}; fixed_thr = 0.95, lost_thr = 0.05, keep_unfinished=false) where A
	X,Y,N,αβ = frequency_series(posh)
	out = Array{FrequencyTraj{A}, 1}(undef, 0)
	for a in 1:size(Y,2)
		#
		start_indices = findall(collect(2:length(X)-1)) do i 
			if Y[i-1,a] >= fixed_thr && Y[i,a] >= fixed_thr # Was fixed
				if Y[i+1,a] < fixed_thr # Now not fixed
					return true
				else
					return false
				end
			elseif Y[i-1,a] <= lost_thr && Y[i,a] <= lost_thr # Was lost
				if Y[i+1,a] > lost_thr  # Now exists
					return true
				else
					return false
				end
			else
				return false
			end
		end .+ 1
		#
		end_indices = findall(collect(2:length(X)-1)) do i 
			if Y[i,a] >= fixed_thr && Y[i+1,a] >= fixed_thr # Now fixed
				if Y[i-1,a] < fixed_thr # Was not fixed
					return true
				else
					return false
				end
			elseif Y[i,a] <= lost_thr && Y[i+1,a] <= lost_thr # Now lost
				if Y[i-1,a] > lost_thr # Was not lost
					return true
				else
					return false
				end
			else
				return false
			end
		end .+ 1


		# Fixing problems with boundaries : trajectories that we did not see start 
		while !isempty(end_indices) && !isempty(start_indices) && end_indices[1] < start_indices[1] 
			splice!(end_indices,1)
		end
		# ... or that we will not see end.
		if !isempty(end_indices) && !isempty(start_indices) && end_indices[end] < start_indices[end]
			if keep_unfinished
				push!(end_indices, length(X)-1)
			elseif !isempty(end_indices) && !isempty(start_indices)
				splice!(start_indices, length(start_indices))
			end
		elseif isempty(end_indices) && !isempty(start_indices)
			if length(start_indices)>1
				error("Empty `end_indices` and more than one in `start_indices`")
			elseif keep_unfinished
				push!(end_indices, length(X)-1)
			else
				start_indices=Int64[]
			end
		end

		# If nothing found, break
		if isempty(start_indices) || isempty(end_indices)
			continue
		end
		# Storing found trajectory
		for (is,ie) in zip(start_indices, end_indices)
			val = αβ[a]
			date = X[is]
			t = X[is:ie] .- X[is]
			freq = Y[is:ie,a]
			pop = N[is:ie]
			strains = [Array{AbstractString, 1}(undef, x) for x in pop]
			index = Dict(:start=>1, :end=>ie-is+1	, :active=>missing)
			fixation = :poly
			if Y[ie,a] >= fixed_thr
				fixation = :fixed
			elseif Y[ie,a] <= lost_thr
				fixation = :lost
			end
			push!(out, FrequencyTraj(posh.i, val, date, t, freq, pop, strains, index, fixation, Dict()))
		end
	end
	return out
end
"""
	all_trajectories(posh::Array{<:PosEvo,1}; fixed_thr = 0.95, lost_thr = 0.05, keep_unfinished=false)
"""
function all_trajectories(posh::Array{<:PosEvo,1}; fixed_thr = 0.95, lost_thr = 0.05, keep_unfinished=false)
	return cat(all_trajectories.(posh, fixed_thr=fixed_thr, lost_thr=lost_thr, keep_unfinished=keep_unfinished)..., dims=1)
end

"""
	previous_state_condition(traj::Array{<:FrequencyTraj,1}, past_state::Symbol)
"""
function previous_state_condition(traj::Array{<:FrequencyTraj,1}, past_state::Symbol)
	idx = zeros(Int64, 0)
	for (id, ft) in enumerate(traj)
		if (ft.freq[1] < 0.5 && past_state==:lost) || (ft.freq[1] > 0.5 && past_state == :fixed)
			push!(idx,id)
		end
	end
	return traj[idx]
end

"""
	frequency_condition(traj::Array{FrequencyTraj{A},1}, α; dα=0.05, shift_time=true) where A
"""
function frequency_condition(traj::Array{FrequencyTraj{A},1}, α; dα=0.05, shift_time=true) where A
	traj_ = Array{FrequencyTraj{A},1}(undef,0)
	for (id, ft) in enumerate(traj)
		for (i,(t,f)) in enumerate(zip(ft.t, ft.freq))
			if f >= α - dα && f < α + dα
				ft_ = deepcopy(ft)
				if shift_time
					ft_.t = ft_.t .- t
					ft_.date = ft_.date + t
					ft_.index[:active] = i
				end
				push!(traj_, ft_)
				break
			end
		end
	end
	return traj_
end
"""
	min_frequency_condition(traj::Array{<:FrequencyTraj,1}, α; shift_time=true)

Return trajectories that reach at least frequency `α`. 
"""
function min_frequency_condition(traj::Array{<:FrequencyTraj,1}, α; shift_time=true)
	dα = (1-α)/2
	return frequency_condition(traj, α+dα, dα=dα, shift_time=shift_time)
end

"""
	population_size_condition(traj::Array{FrequencyTraj,1}, minpop)

Select trajectories that are based on a total population of at least `minpop` for their whole duration.  
Values of `mode`: `:overall`, `active` or `:end`.
"""
function population_size_condition(traj::Array{<:FrequencyTraj,1}, minpop; mode=:overall)
	idx = zeros(Int64, 0)
	for (id, ft) in enumerate(traj)
		flag = true
		if mode == :overall
			for p in ft.pop[2:end-1]
				if p < minpop
					flag = false
				end
			end
		elseif mode == :active
			flag = ft.pop[ft.index[:active]] >= minpop
		elseif mode == :end
			flag = ft.pop[end] >= minpop
		else
			@error "Unrecognized mode $mode"
		end
		if flag
			push!(idx, id)
		end
	end
	return traj[idx]
end

"""
	derivative_condition(traj::Array{<:FrequencyTraj,1}; direction=:positive)

Note: Use the `:active` index of trajectories for computing derivative. Zero derivative are ignored. 
"""
function derivative_condition(traj::Array{<:FrequencyTraj,1}; direction=:positive)
	idx = Int64[]
	for (id, ft) in enumerate(traj)
		df = (ft.index[:active] < 2) ? 0. : (ft.freq[ft.index[:active]] - ft.freq[ft.index[:active] - 1])
		if direction == :positive && df > 0
			push!(idx, id)
		elseif direction == :negative && df < 0
			push!(idx, id)
		end
	end
	return traj[idx]
end


"""
	get_strains!(traj::FrequencyTraj, fp::FluPop)
	get_strains!(traj::Array{<:FrequencyTraj}, fp::FluPop)

Find strains corresponding to trajectory `traj`, filling the `traj.strains` field.
"""
function get_strains!(traj::FrequencyTraj, fp::FluPop)
	for (i,t) in enumerate(traj.t)
		db = find_datebin(traj.date + t, fp)
		traj.strains[i] = [x[:strain] for x in find_strains(fp.datebin[db], traj.i, traj.val)]
	end
end
function get_strains!(traj::Array{<:FrequencyTraj}, fp::FluPop)
	for t in traj
		get_strains!(t, fp)
	end
end

"""
	compute_tree_spread!(traj::FrequencyTraj, tree::Tree)

Find how `traj` is spread in the tree: does come from one clade, two, more? This is quantified using either entropy or the sum sum of squared frequencies for the distribution of `traj` into clades. 

## Notes
- This only makes sense for mutations that are *rising*. Mutations falling from a fixed position to a finite frequency will not appear in the tree at the moment they start falling! 
- This function assumes `traj.strains` is already known. If it hasn't been determined before, call compute_tree_spread!(traj::FrequencyTraj, tree::Tree, fp::FluPop)
"""
function compute_tree_spread!(traj::FrequencyTraj, tree::Tree)
	traj.data[:treespread] = [Int64[] for i in 1:length(traj)]
	for (i, strains) in enumerate(traj.strains)
		class = Dict()
		for label in strains
			if haskey(tree.lleaves, label)
				# Find branch upstream of given strain that bears mutation
				r = TreeTools.find_mut_root(tree.lleaves[label], traj.i, traj.val)
				class[r] = get(class, r, 0) + 1
			end
		end
		traj.data[:treespread][i] = collect(values(class))
	end
end
"""
	compute_tree_spread!(traj::FrequencyTraj, tree::Tree, fp::FluPop)
"""
function compute_tree_spread!(traj::FrequencyTraj, tree::Tree, fp::FluPop)
	get_strains!(traj, fp)
	compute_tree_spread!(traj, tree)
end

"""
	get_regions!(traj::FrequencyTraj, fp::FluPop)
	get_regions!(traj::Array{<:FrequencyTraj}, fp::FluPop) 
"""
function get_regions!(traj::FrequencyTraj, fp::FluPop)
	get_strains!(traj, fp)
	traj.data[:regions] = [Dict{String,Int64}() for i in 1:length(traj)]
	traj.data[:countries] = [Dict{String,Int64}() for i in 1:length(traj)]
	for (i,strains) in enumerate(traj.strains)
		# traj.regionspread[i] = Dict(r=>0 for r in nextstrain_regions)
		for s in strains
			!in(fp.strains[s][:region], ["?",'?']) && (traj.data[:regions][i][fp.strains[s][:region]] = get(traj.data[:regions][i], fp.strains[s][:region], 0) + 1 )
			!in(fp.strains[s][:country], ["?",'?']) && (traj.data[:countries][i][fp.strains[s][:country]] = get(traj.data[:countries][i], fp.strains[s][:country], 0) + 1 )
		end
	end
end
function get_regions!(traj::Array{<:FrequencyTraj}, fp::FluPop)
	for t in traj
		get_regions!(t, fp)
	end
end

"""
	get_date_idx(traj::FrequencyTraj, date::Date)
	get_date_idx(traj::FrequencyTraj, date)

Index at which `date` appears in `traj`. Return `nothing` if nothing is found. `date` should be convertable to a `Date` type.
"""
function get_date_idx(traj::FrequencyTraj, date::Date)
	for i in 1:length(traj)
		if traj.date + traj.t[i] == date
			return i
		end
	end
	return nothing
end
get_date_idx(traj::FrequencyTraj, date) = get_date_idx(traj, Date(date))

"""
	find_common_dates(t1::FrequencyTraj, t2::FrequencyTraj)

Return two vectors of indices corresponding to the common dates of `t1` and `t2`. 
"""
function find_common_dates(t1::FrequencyTraj, t2::FrequencyTraj)
	if t1.date > t2.date 
		idx2, idx1 = find_common_dates(t2, t1)
		return idx1, idx2
	end
	idx1 = Int64[]; idx2 = Int64[]
	for i1 in 1:length(t1)
		i2 = get_date_idx(t2, t1.date + t1.t[i1])
		if !isnothing(i2) 
			push!(idx1, i1)
			push!(idx2, i2)
		end
	end
	return idx1, idx2
end

"""
	traj_distance(t1::FrequencyTraj, t2::FrequencyTraj; dist_type=:l1, kwargs...)

Distance types are `(:l1, :linfit)`. 
"""
function traj_distance(t1::FrequencyTraj, t2::FrequencyTraj; dist_type=:strains, kwargs...)
	!in(dist_type, traj_distance_types()) && @error("Unknown distance type $dist_type.")
	if dist_type == :l1
		return traj_distance_l1(t1, t2; kwargs...)
	elseif dist_type == :linfit
		return traj_distance_linfit(t1, t2; kwargs...)
	elseif dist_type == :strains
		return traj_distance_strains(t1,t2; kwargs...)
	end
end

function traj_distance_default_kwargs(dist_type::Symbol)
	!in(dist_type, traj_distance_types()) && @error("Unknown distance type $dist_type.")
	if dist_type == :l1
		return (normalize=:full_length, common_only=false)
	elseif dist_type == :linfit
		return (min_overlap=3,)
	elseif dist_type == :strains
		return (normalize=true, min_overlap=3)
	end
end
traj_distance_types() = (:l1, :linfit, :strains)

"""
	traj_distance_strains(t1::FrequencyTraj, t2::FrequencyTraj; min_overlap=3)
"""
function traj_distance_strains(t1::FrequencyTraj, t2::FrequencyTraj; normalize=true, min_overlap=3)
	if t1.date > t2.date
		return traj_distance_strains(t2, t1, min_overlap=min_overlap, normalize=normalize)
	end
	# 
	idx1, idx2 = find_common_dates(t1,t2)
	if length(idx1) < min_overlap 
		return 1.
	end	
	# 
	dif = 0.
	Z = 0.
	for (i1,i2) in zip(idx1, idx2)
		# if isempty(t1.strains[i1]) || isempty(t2.strains[i2])
		# 	@error("Trajectory without strain :$(t1.i)$(t1.val)")
		# end
		dif += 1 - jaccard(t1.strains[i1], t2.strains[i2])
		if !isempty(t1.strains[i1],) || !isempty(t2.strains[i2],)
			Z += 1
		end
	end
	if normalize
		return dif / Z
	else
		return dif
	end
end
function jaccard(A,B) 
	try
		if isempty(A) && isempty(B)
			return 1.
		elseif !isempty(A) && !isempty(B)
			length(intersect(A,B)) / length(union(A,B))
		else
			return 0.
		end
	catch err
		println(A)
		println(B)
		@error(err)
	end
end

"""
	traj_distance_linfit(t1::FrequencyTraj, t2::FrequencyTraj; min_overlap=3)
"""
function traj_distance_linfit(t1::FrequencyTraj, t2::FrequencyTraj; min_overlap=3)
	if t1.date > t2.date
		return traj_distance_linfit(t2,t1,min_overlap=min_overlap)
	end
	# 
	idx1, idx2 = find_common_dates(t1,t2)
	if length(idx1) < min_overlap 
		return 1.
	end
	# 
	f1 = t1.freq[idx1]
	f2 = t2.freq[idx2]
	line = Polynomials.fit(f1,f2,1)
	return round(abs(line[1] - 1.), digits=5)
end
"""
	traj_distance_l1(t1::FrequencyTraj, t2::FrequencyTraj; normalize=true, common_only=false)

Sum `abs(f1-f2)` for all the length of trajectories `t1` and `t2`. When only one trajectory is active, the value taken for the other is either 0. or 1. 
"""
function traj_distance_l1(t1::FrequencyTraj, t2::FrequencyTraj; normalize=false, common_only=false)
	if t1.date > t2.date
		return traj_distance_l1(t2,t1,normalize=normalize, common_only=common_only)
	end
	# 
	idx1, idx2 = find_common_dates(t1,t2)
	# If common support is null
	α = 0.
	if length(idx1) <= α * length(t1) || length(idx2) <= α * length(t2)
		if normalize 
			return 1
		else
			return max(length(t1), length(t2))
		end
	end
	#
	f2 = t2.freq[1] < 0.5 ? 0. : 1. # Initial value of t2
	dif = 0.
	Z = 0 # Normalization
	# Before t2 starts
	if !common_only
		i1 = 1
		while i1 < idx1[1]
			dif += abs(t1.freq[i1] - f2)
			Z += 1
			i1 += 1
		end
	end

	# Common part of both trajectories
	for (i1,i2) in zip(idx1, idx2)
		dif += abs(t1.freq[i1] - t2.freq[i2])
		Z += 1
	end

	# Continue counting differences if one of the trajectories has ended and not the other
	if !common_only
		i1 = idx1[end]; i2 = idx2[end]
		while i1 < length(t1)-1 || i2 < length(t2)-1
			i1 += 1
			i2 += 1
			if i1 > length(t1)
				f1 = t1.fixation==:fixed ? 1. : 0.
			else
				f1 = t1.freq[i1]
			end
			if i2 > length(t2)
				f2 = t2.fixation==:fixed ? 1. : 0.
			else
				f2 = t2.freq[i2]
			end
			dif += abs(f1 - f2)
			Z += 1
		end
	end
	# 
	if !normalize
		Z = 1
	end
	return sum(diff) / Z
end

"""
	compute_distance_matrix(f, trajectories)

Compute pairwise distance matrix of `trajectories` for distance function `f`. 
"""
function compute_distance_matrix(f, trajectories)
	dist_mat = zeros(Float64, length(trajectories), length(trajectories))
	for (i,t1) in enumerate(trajectories)
	    for j in (i+1):length(trajectories) 
	        dist_mat[i,j] = f(t1,trajectories[j])	
	        dist_mat[j,i] = dist_mat[i,j]
	    end
	    dist_mat[i,i] = 0.
	end	
	return dist_mat
end

"""
	cluster_trajectories_naive(f, trajectories; thr=0.)

Group together trajectories that are exactly the same. Distance function is `f`. 
"""
function cluster_trajectories_naive(f, trajectories; thr=0.)
	dist_mat = compute_distance_matrix(f, trajectories)
	return cluster_trajectories_naive(dist_mat, trajectories, thr=thr)
end
"""
	cluster_trajectories_naive(dist_mat::Array{<:Real,2}, trajectories; thr=0.)
"""
function cluster_trajectories_naive(dist_mat::Array{<:Real,2}, trajectories; thr=0.)
	clusters = Dict{Int64,Array{Int64,1}}()
	clustered = zeros(Bool, length(trajectories)) # Is the trajectory in a cluster 
	for i in 1:length(trajectories)
	    clustered[i] && continue
	    clusters[i] = Int64[i]
	    clustered[i] = true
	    increase_cluster_naive!(clusters[i], clustered, dist_mat, thr)
	end
	return clusters
end
function increase_cluster_naive!(cluster, clustered::Array{Bool,1}, dist_mat, thr::Real)
	cflag = true
	idx = 1
	while cflag
		cflag = false
		for i in idx:length(cluster)
			for j in 1:size(dist_mat,1)
				if !clustered[j] && dist_mat[cluster[i],j] <= thr # If j is not clustered already and is close to `cluster[i]`
					push!(cluster, j)
					clustered[j] = true
					cflag = true
				end
			end
		end
		idx += 1
	end
end

function traj_clustering_default_kwargs(clustering_type::Symbol)
	if clustering_type == :naive
		return (thr=0.,)
	else
		@error("Unknown clustering type $clustering_type")
	end
end

"""
	effective_trajectories(trajectories; clustering=:naive, distance=:l1, common_only=true)

Cluster trajectories based on similarity, and return one representative per cluster
"""
function effective_trajectories(trajectories; 
	clustering=:naive, clust_kwargs=traj_clustering_default_kwargs(clustering),
	distance=:strains, dist_kwargs=traj_distance_default_kwargs(distance))

	# Could vary based on kwargs
	dfunc(t1, t2) = traj_distance(t1,t2; dist_type=distance, dist_kwargs...)
	# clusters -- could also vary based on kwargs
	clusters = cluster_trajectories_naive(dfunc, trajectories; clust_kwargs...)

	# Picking one representative per cluster
	efftraj = []
	for (i,idx) in clusters
		push!(efftraj, trajectories[idx[rand(1:length(idx))]])
	end
	return efftraj
end
