export alphabet

"""
	abstract type AbstractStrain end

The only requirement for concrete implementations of `AbstractStrain` is the existence of two fields: `seq` and `data::Dict`. `data` is often expected to contain a key `"strain"` referring to the name of the strain by default, but I should try not to enforce it. 
"""
abstract type AbstractStrain end

getindex(s::AbstractStrain, key::String) = get(s.data, key) do 
	haskey(s.data, Symbol(key)) ? s.data[Symbol(key)] : @error "Keys $key or $(Symbol(key)) not found."
end
getindex(s::AbstractStrain, key::Symbol) = get(s.data, key) do 
	haskey(s.data, String(key)) ? s.data[String(key)] : @error "Keys $key or $(String(key)) not found."
end
"""
	isempty(st::AbstractStrain)
"""
function isempty(st::AbstractStrain)
	return isempty(st.seq) && isempty(st.data)
end

"""
	mutable struct Strain{A} <: AbstractStrain

Fields: 
- `seq::BioSequence{A}`
- `data::Dict`
"""
mutable struct Strain{A} <: AbstractStrain
	seq::LongSequence{A}
	data::Dict
end
function Strain(seqtype)
	if seqtype == :aa
		return Strain(LongAminoAcidSeq(), Dict())
	elseif seqtype == :dna
		Strain(LongDNASeq(), Dict())
	elseif seqtype == :rna
		Strain(LongRNASeq(), Dict())
	else
		unknown_seqtype()
	end
end
"""
	Strain(seq, dat, seqtype)

Construct `Strain` object from `seq` and `dat`, when `seq` needs to be converted to a `BioSequence` object first. Values of `seqtype`: `aa`, `dna` or `rna`. 
"""
function Strain(seq, dat, seqtype ; verbose=false)
	if seqtype == :aa
		try 
			return Strain(LongAminoAcidSeq(seq), dat)
		catch
			cannot_read(seq, dat, seqtype)
			return Strain(seqtype)
		end
	elseif seqtype == :dna
		try
			return Strain(LongDNASeq(seq), dat)
		catch
			cannot_read(seq, dat, seqtype)
			return Strain(seqtype)
		end
	elseif seqtype == :rna
		return Strain(LongRNASeq(seq), dat)
	else
		unknown_seqtype()
	end
end
function cannot_read(seq, dat, seqtype)
	@warn "Could not read the following sequence:"
	println("$seq")
	println("Associated data:\n$dat")
	println("Called seq. type: $seqtype")
	println()
	return nothing
end

"""
	mutable struct ArtificialStrain{T} <: AbstractStrain
		seq::Array{T,1}
		data::Dict
	end
"""
mutable struct ArtificialStrain{T} <: AbstractStrain
	seq::Array{T,1}
	data::Dict
end
function ArtificialStrain(T::DataType)
	return ArtificialStrain(T[], Dict())
end
"""
	ArtificialStrain(seq, dat, T::DataType = Int64)
"""
function ArtificialStrain(seq, dat, T::DataType = Int64)
	return ArtificialStrain(convert(Array{T,1}, seq), dat)
end


# abstract type SequenceProfile end
# For now I don't see an efficient application of profile outside of this. 
# For DCA like stuff, it's a float vector, and it's highly unpractical to have something both general and efficient here. 
# So no abstract type

# One way would be to have `data` be an array of subtype `colProfile`.
# `colProfile` could then be designed independently
# As long as the getindex function is overloaded properly this could be practical.
"""
	mutable struct SiteFrequency{A}
		i::Int64
		M::Int64
		freq::Dict{A,Float64}
	end

## Note
- From outside, the type should behave like a `Dict{A,Float64}`. That is, `keys(F::SiteFrequency)` should return the keys of `F.freq`. 
- TRY: Calling `f[a]` returns `0.` as a default value if `a` does not appear in `f.freq`. 
"""
mutable struct SiteFrequency{A}
	i::Int64
	M::Int64
	freq::Dict{A,Float64}
end
"""
	SiteFrequency(T::DataType, i::Int64=0, M::Int64=0)

Constructor for empty structure. 
"""
function SiteFrequency(T::DataType, i::Int64=0, M::Int64=0)
	return SiteFrequency(i, M, Dict{T,Float64}())
end
setindex!(F::SiteFrequency{A}, f::Float64, a::A) where A = setindex!(F.freq, f, a)
get(F::SiteFrequency{A}, a::A, def) where A = get(F.freq, a, def)
getindex(F::SiteFrequency{A}, a::A) where A = get(F, a, 0.)
count(F::SiteFrequency{A}, a::A) where A = (F.M == 0) ? error("SiteFrequency: Unknown sequence number") : round(Int64, F.freq[a] * F.M)
counts(F::SiteFrequency) = (F.M == 0) ? error("SiteFrequency: Unknown sequence number") : Dict(k=>v*F.M for (k,v) in F.freq)
keys(F::SiteFrequency) = keys(F.freq)
alphabet(F::SiteFrequency) = keys(F)
values(F::SiteFrequency) = values(F.freq)
isempty(F::SiteFrequency) = (F.M == 0) || isempty(F.freq)

"""
	mutable struct Profile{A}
		data::Array{Dict{A, Float64},1}
		M::Int64
	end
Access `data` by indexing: `p[i,a]` --> `p.data[i][a]`, or `p[i][a]` --> `p.data[i][a]`. This does not return a default value. For a default value, try `get(p, i, a, def)` or `get(p[i], a, def)` (the former might be slower). 
"""
mutable struct Profile{A}
	data::Array{SiteFrequency{A},1}
	M::Int64 # Number of sequences the profile is based on
end
function Profile(T::DataType, L::Int64, M::Int64=0)
	return Profile([SiteFrequency(T,i,M) for i in 1:L], M)
end
getindex(P::Profile, i::Int64, a) = P.data[i][a]
getindex(P::Profile, i::Int64) = P.data[i]
get(P::Profile, i::Int64, default) = get(P.data, i, default)
get(P::Profile, i::Int64, a, default) = get(get(P.data, i, Dict()), a, default)
length(P::Profile) = length(P.data)
iterate(P::Profile, n=1) = iterate(P.data, n)

"""
	Profile(S::Array{<:BioSequence{A}}; ambiguous=false) where A
	Profile(S::Array{<:Strain{A}}; ambiguous=false) where A
"""
function Profile(S::Array{<:BioSequence{A}}; ambiguous=false) where A
	if isempty(S)
		@error "Cannot build a profile  from empty alignment"
	else
		prof = Profile(eltype(A), length(S[1]), length(S))
		for (m,s) in enumerate(S)
			for (i,a) in enumerate(s)
				if ambiguous || !isambiguous(a)
					prof.data[i][a] = get(prof[i], a, 0.) + 1. 
				else
					prof.data[i].M -= 1
				end
			end
		end
		# Normalizing 
		for i in 1:length(prof)
			for a in alphabet(prof.data[i])
				prof.data[i][a] /= prof.data[i].M
			end
		end
	end
	return prof
end
function Profile(S::Array{<:Strain{A}}; ambiguous = false) where A
	M = length(S)
	if isempty(S)
		@error "Cannot build a profile from empty alignment"
	else
		prof = Profile(eltype(A), length(S[1].seq), length(S))
		for (m,s) in enumerate(S)
			for (i,a) in enumerate(s.seq)
				if ambiguous || !isambiguous(a)
					prof.data[i][a] = get(prof[i], a, 0.) + 1. 
				else
					prof.data[i].M -= 1
				end
			end
		end
		# Normalizing
		for i in 1:length(prof)
			for a in alphabet(prof.data[i])
				prof.data[i][a] /= prof.data[i].M
			end
		end
	end
	return prof
end
