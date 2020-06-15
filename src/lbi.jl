export get_lbi_timescale, get_lbi!

"""
	get_lbi_timescale(lineage, segment)

Return the timescale used for computing LBI values, in branch length units.
 """
function get_lbi_timescale(lineage, segment)
	if typeof(lbi_integration_time) == Year
		ct = lbi_integration_time.value 
	elseif typeof(lbi_integration_time) == Month
		ct = lbi_integration_time.value / 12
	elseif typeof(lbi_integration_time) == Day
		ct = lbi_integration_time.value / 365
	else
		ct = lbi_integration_time
	end
	τ = substitution_rate[lineage, segment] * ct
	return τ
end


"""
	get_lbi!(t::Tree{LBIData}, fp::FluPop, strainnames::Array{String}, datemin, datemax, τ, lbi_field)

Compute lbi of `strainsnames` according to tree `t`. Leaves of tree posterior to `datemax` or anterior to `datemin` are set to dead. Tree is modified in the process (life-state and lbi value of nodes). 
Store the resulting lbi value in `fp.strains[lbi_field]`
"""
function get_lbi!(t::Tree{LBIData}, fp::FluPop, strainnames::Array{<:AbstractString}, datemin, datemax, τ, lbi_field, verbose=false)
	# Setting life-state of tree nodes
	for node in values(t.leaves)
		if haskey(fp.strains, node.label)
			node.data.alive = datemin <= fp.strains[node.label][:date] < datemax 
		else
			node.data.alive = false
		end
	end
	TreeTools.set_live_nodes!(t, set_leaves=false)	
	# Computing LBI
	TreeTools.lbi!(t, τ, normalize=true) && verbose && begin
		println("LBI computation failed")
	    println("-- ", length(strainnames), " strains")
	    println("-- Date range $datemin - $datemax")
	end
	for (i,s) in enumerate(strainnames)
		if haskey(t.lleaves, s)
			fp.strains[s].data[lbi_field] = t.lleaves[s].data.lbi
		else
			fp.strains[s].data[lbi_field] = missing 
		end
	end
	nothing
end

"""
	function get_lbi!(fp::FluPop, t::Tree{LBIData}, datestyle=:all_anterior; 
				lineage="h3n2", 
				segment="ha", 
				τ=get_tau(lineage,segment),
				lbi_field=:lbi)

Compute lbi of strains in `fp` based on tree `t` and date binning of `fp.datebin`. Tree is modified in the process (life-state of nodes).
The `datestyle` arguments modifies the life-state of nodes when LBI is computed. If `:all_anterior`, all nodes anterior to the current datebin are set to live. If `:datebin`, only nodes in the current datebin are set to live. 
If a strain does not exist in the tree, its lbi will be set to `missing`. 

"""
function get_lbi!(fp::FluPop, t::Tree{LBIData}, datestyle=:all_anterior; 
				lineage="h3n2", 
				segment="ha", 
				τ=get_lbi_timescale(lineage,segment),
				lbi_field=:lbi, 
				verbose=false)
	# Find `datemin` as a function of `datestyle` 
	if datestyle == :all_anterior
		datemin = let val = findmin(collect(keys(fp.datebin)))[1][1]
			(x) -> val
		end
	elseif datestyle == :datebin
		datemin = (x) -> minimum(x)
	else
		error("Accepted values for `datestyle`: `:all_anterior` or `:datebin`")
	end
	# Compute LBI for all datebins
	for (datebin,strains) in sort(OrderedDict(fp.datebin), by=x->x[1])
		if length(strains)>1
			labels = [x[:strain] for x in strains]
			get_lbi!(t, fp, labels, datemin(datebin), max(datebin...), τ, lbi_field, verbose)
		else
			for s in strains
				s.data[lbi_field] = missing
			end
		end
	end
	# Setting lbi to missing for strains not in a datebin
	for s in values(fp.strains)
		if !haskey(s.data, lbi_field)
			s.data[lbi_field] = missing
		end
	end
end


"""
	local_lbi_maximas(t::Tree{LBIData}, fp::FluPop, τ)

The idea is to find local maximas of the LBI for each date bin of `fp`. In order to determine if `n` is a local maximum, one needs to know the LBI of its ancestor `a`. Thus, the lbi of `a` has to be recomputed for every datebin, when new leaves come live and old ones die. 
## Notes
- `fp` should be date-binned already. 
- Strains in `fp` and `t` should match **exactly**. 
"""
function local_lbi_maximas!(t::Tree{LBIData}, fp::FluPop, τ=get_lbi_timescale("h3n2","ha"))
	out = Dict()
	# Setting nodes alive in date order
	for (db, strains) in sort(fp.datebin, by=x->x[1])
		for n in values(t.lnodes)
			n.data.alive = false
		end
		for s in strains
			t.lleaves[s[:strain]].data.alive = true
		end
		TreeTools.set_live_nodes!(t, set_leaves=false)	
		TreeTools.lbi!(t, τ, normalize=true)
		# Local maximas for the current tree
		llm = local_lbi_maximas(t)
		out[db] = [(label=llm[i], lbi=t.lnodes[llm[i]].data.lbi) for i in 1:length(llm)]
	end
	return out
end
"""
	local_lbi_maximas(t::Tree{LBIData})

Find all nodes that are a local maxima of the LBI. Does not compute the LBI. 
"""
function local_lbi_maximas(t::Tree{LBIData})
	out = String[]
	for (l,n) in t.lnodes
		if n.isroot || n.data.lbi > n.anc.data.lbi 
			if mapreduce(x-> n.data.lbi > x.data.lbi, *, n.child, init=true)
				push!(out, n.label)
			end
		end
	end
	return out
end


# if false 

# """
# 	local_lbi_maximas(t::Tree{LBIData}, sp::StrainPop, τ)

# The idea is to find local maximas of the LBI for each date bin of `sp`. In order to determine if `n` is a local maximum, one needs to know the LBI of its ancestor `a`. Thus, the lbi of `a` has to be recomputed for every datebin, when new leaves come live and old ones die. 
# ## Notes
# - `sp` should be date-binned already. 
# - Strains in `sp` and `t` should match **exactly**. 
# """
# function local_lbi_maximas!(t::Tree{LBIData}, sp::StrainPop, τ)
# 	out = Dict()
# 	# Setting nodes alive in date order
# 	for (db, strains) in sort(sp.datebin, by=x->x[1])
# 		for n in values(t.lnodes)
# 			n.data.alive = false
# 		end
# 		for s in strains
# 			t.lleaves[s.strain].data.alive = true
# 		end
# 		TreeTools.set_live_nodes!(t, set_leaves=false)	
# 		TreeTools.lbi!(t, τ, normalize=true)
# 		# Local maximas for the current tree
# 		llm = local_lbi_maximas(t)
# 		out[db] = [(label=llm[i], lbi=t.lnodes[llm[i]].data.lbi) for i in 1:length(llm)]
# 	end
# 	return out
# end
# """
# 	local_lbi_maximas(t::Tree{LBIData})

# Find all nodes that are a local maxima of the LBI. Does not recompute LBI. 
# """
# function local_lbi_maximas(t::Tree{LBIData})
# 	out = String[]
# 	for (l,n) in t.lnodes
# 		if n.isroot || n.data.lbi > n.anc.data.lbi 
# 			if mapreduce(x-> n.data.lbi > x.data.lbi, *, n.child, init=true)
# 				push!(out, n.label)
# 			end
# 		end
# 	end
# 	return out
# end

# end