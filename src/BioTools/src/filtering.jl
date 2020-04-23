"""
	remove(strains::Array{<:AbstractStrain}, names, [field=:strain]; verbose=false)

Remove strains for which `s.data[field]` is in `names`. 
"""
function remove(strains::Array{<:AbstractStrain}, names, field=:strain; verbose=false)
	out = Array{eltype(strains),1}(undef, length(strains))
	idx = 1
	for s in strains
		if !in(s[field], names)
			out[idx] = s
			idx += 1
		end 
	end
	verbose && println("Filtered $(length(strains) - idx + 1) strains")
	deleteat!(out, idx:length(strains))
end
"""
	remove!(strains::Array{<:AbstractStrain}, names, [field=:strain]; verbose=false)

Remove strains for which `s.data[field]` is in `names`. 
"""
function remove!(strains::Array{<:AbstractStrain}, names, field=:strain; verbose=false)
	indices = Int64[]
	for (i,s) in enumerate(strains)
		if in(s[field], names)
			push!(indices, i)
		end 
	end	
	verbose && println("Filtered $(length(indices)) strains")
	deleteat!(strains, indices)
end

"""
	gapfilter(s::Strain; threshold=0.1)

Return true if `s.seq` has less than `threshold * length(s.seq)` gaps.
"""
gapfilter(s::Strain; threshold=0.05) = countgaps(s.seq) < threshold * length(s.seq)
"""
	gapfilter(s::ArtificialStrain; threshold=0.)
"""
gapfilter(s::ArtificialStrain; threshold=0.)  = true
# gapfilter(threshold::Float64) = s->gapfilter(s, threshold)

hasdate(x::AbstractStrain) = (!ismissing(get(x.data, "date", missing)) || !ismissing(get(x.data, :date, missing)))
