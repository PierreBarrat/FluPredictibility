"""
	FluPop(f::Union{AbstractString,IO}, sequencetype::Symbol, headerfields; 
	flulineage=missing, segment=missing, strainfilters = [!is_flu_outlier(flulineage), BioTools.hasdate], separator = '|')

Call `readfastastrains` to read `f`. Store the result in a `FluPop` object. 
"""
function FluPop(f::Union{AbstractString,IO}, 
	sequencetype, 
	headerfields; 
	flulineage=missing, 
	segment=missing, 
	strainfilters = [!is_flu_outlier(flulineage), BioTools.hasdate, s->BioTools.gapfilter(s,threshold=0.05)],
	separator = '|', 
	ignore_read_errors=false)	
	
	strains = readfastastrains(f, sequencetype, headerfields, separator = separator, strainfilters = strainfilters, ignore_read_errors=ignore_read_errors)
	return FluPop(strains = Dict{String, eltype(strains)}(String(x[:strain])=>x for x in strains))
end
"""
	AAFluPop(f::Union{AbstractString,IO}, headerfields; [kwargs...])

Call `FluPop(f, :aa, headerfields, [kwargs...])
"""
function AAFluPop(f::Union{AbstractString,IO}, headerfields; 
	flulineage=missing, 
	segment=missing, 
	strainfilters = [!is_flu_outlier(flulineage), BioTools.hasdate],
	separator = '|')
	return FluPop(f, :aa, headerfields, flulineage=flulineage, segment=segment, strainfilters=strainfilters, separator=separator)
end

is_flu_outlier(x::AbstractStrain, flulineage) = in(x[:strain], outliers[flulineage])
function is_flu_outlier(flulineage)
	return x->is_flu_outlier(x,flulineage)
end

"""
	read_mutations!(tree::Tree, mutfile::String)

Read mutations from a JSON file outputted by augur. Store them into `tree`. 

## Note
This is in the `Flu` module because it is explicitely designed for use of augur on flu. In particular, it uses the `gene_positions` global. 
"""
function read_mutations!(tree::Tree, mutfile::String; type=:aa, lineage="h3n2", segment="ha")
	if type == :aa
		read_mutations_aa!(tree, mutfile, lineage, segment)
	else
		@warn "Only implemented for `type=:aa`"
	end
end

"""
"""
function read_mutations_aa!(tree, mutfile::String, lineage, segment)
	muts = JSON.Parser.parsefile(mutfile)["nodes"]
	for (label, mut) in muts
		if haskey(tree.lnodes, label)
			tmp = Array{TreeTools.Mutation}(undef, 0)
			for (gene, pos) in Flu.gene_positions[lineage, segment]
				if !haskey(mut["aa_muts"], "HA1")
					println(label); 
				else
					for m in mut["aa_muts"][gene]
						push!(tmp, _parse_aa_mut(m))
						tmp[end].i = tmp[end].i + Int64((pos[1] - 1)/3) # Offset for different genes
					end
				end
			end
			tree.lnodes[label].data.mutations = tmp
		end
	end
end
"""
Parse a string of format `XiY` into a `TreeTools.Mutation` object with fields `i`, `X` and `Y`.
"""
function _parse_aa_mut(m::String)
	if length(m) < 3
		error("Can't parse mutation string of length smaller than 3")
	end
	i = parse(Int64, m[2:end-1])
	old = AminoAcid(m[1])
	new = AminoAcid(m[end])
	out =  TreeTools.Mutation(i, old, new)
	return out
end
