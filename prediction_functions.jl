using TreeTools
using JSON

"""
Split strains in clades using a timetree. Two leaves `n1` and `n2` belong to the same clade if and only if `lca(n1,n2)` lives after `datethr`. 
"""
function split_in_clades(t::Tree{LBIData}, strainnames::Array{String,1}, datethr)
	clades = Dict()
	for s in strainnames
		found = false
		for id in keys(clades)
			if lca(t.lleaves[s], t.lleaves[id]).data.date >= datethr
				push!(clades[id], s)
				found = true 
				break
			end
		end
		if !found
			clades[s] = [s]
		end
	end
	return clades
end

function split_in_clades(t::Tree{LBIData}, strainnames)
	println("Splitting in clades using Kmeans")
end

"""
Stores the dates contained in parsed json file `json` into `t`. The json file should be the output of `augur refine --timetree` 
"""
function dates_to_tree!(t::Tree{LBIData}, json::Dict)
	for (n,dat) in json["nodes"]
		if haskey(t.lnodes, n) 
			t.lnodes[n].data.date = Date(dat["date"])
		end
	end
end
dates_to_tree!(t::Tree{LBIData}, jsonfile::String) = dates_to_tree!(t, JSON.Parser.parsefile(jsonfile))

"""
Honestly, I see nothing super easy but going through all the leaves. 

If I wanted to actually make this efficient, it would need to be some variation of Dijkstra's algorithm: go in a circle around `r` until you encounter a leaf. Actually it's quite easy on a tree : since you can't loop, you don't have to keep a registry of already visited nodes. 
"""
function closest_tree_leaf(r::TreeNode, t::Tree)
	m = Inf 
	out = ""
	for n in values(t.leaves)
		tmp = TreeTools.node_divtime(r, n)
		if tmp < m
			m = tmp
			out = n.label
		end
	end
	return (out, m)
end
closest_tree_leaf(r::String, t::Tree) = closest_tree_leaf(t.lnodes[r], t)

function closest_tree_leaf(r::TreeNode, t::Tree, candidates::Array{<:AbstractString,1})
	m = Inf 
	out = ""
	for s in candidates
		n = t.lnodes[s]
		tmp = TreeTools.node_divtime(r, n)
		if tmp < m
			m = tmp
			out = n.label
		end
	end		
	return (out, m)
end
closest_tree_leaf(r::String, t::Tree, candidates::Array{<:AbstractString,1}) = closest_tree_leaf(t.lnodes[r], t, candidates)

"""
	cluster_lbi_maxima!(t::Tree, fp)

Cluster sequences of each datebin in `fp.datebin` using local maximas of the LBI for the tree of that datebin. 
"""
function cluster_lbi_maxima!(t::Tree, fp)
	LLM = Flu.local_lbi_maximas!(t, fp);
	clusters = Dict()
	lbis = Dict()
	for (db, llm) in LLM
		centers = [x.label for x in llm]
		clusters[db] = Dict(k=>[] for k in centers)
		lbis[db] = Dict(k=>v for (k,v) in llm)
		pop = fp.datebin[db]
		# For each strain in the current pop, find which LBI maximum it is closest to
		for s in pop
			closest_center = closest_tree_leaf(t.lleaves[s[:strain]], t, centers)[1]
			push!(clusters[db][closest_center], s[:strain])
		end
	end
	for (db,cl) in clusters
		for (k,c) in cl
			if isempty(c)
				delete!(clusters[db], k)
				delete!(lbis[db], k)
			end
		end
	end
	return clusters, lbis
end