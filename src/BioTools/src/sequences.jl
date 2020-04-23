export hamming, profile, consensus, numseq

"""
	countgaps(s::BioSequence)
"""
function countgaps(s::BioSequence)
	return BioSequences.count(x->isgap(x)||isambiguous(x), s)
end
countgaps(s::AbstractStrain) = countgaps(s.seq)

function hamming_ignore_ambiguous(a::BioSequence, b::BioSequence)
	h = 0
	for (x,y) in zip(a,b)
		if x != y && !isambiguous(x) && !isambiguous(y)
			h += 1
		end
	end
	return h
end
function hamming_ambiguous_mismatch(a::BioSequence, b::BioSequence)
	h = 0
	for (x,y) in zip(a,b)
		if x != y 
			h += 1
		end
	end
	return h	
end
"""
	hamming(a::BioSequence, b::BioSequence, [normalization=false]; ambiguous=:ignore)
	hamming(a::AbstractStrain, b::AbstractStrain [, normalization=false]; ambiguous=:ignore)

`ambiguous` keyword:
- `:ignore` : ambiguous symbols are considered equal to anything .
- `:mismatch` : ambiguous symbols behave like regular symbols.
"""
function hamming(a::BioSequence, b::BioSequence, 
	normalization=false; 
	ambiguous=:ignore)
	if ambiguous == :ignore && normalization
		return hamming_ignore_ambiguous(a,b)/length(a)
	elseif ambiguous == :ignore && !normalization
		return Float64(hamming_ignore_ambiguous(a,b))
	elseif ambiguous == :mismatch && normalization
		return hamming_ambiguous_mismatch(a,b)/length(a)
	elseif ambiguous == :mismatch && !normalization
		return Float64(hamming_ambiguous_mismatch(a,b))
	else
		@error "Incorrect arguments."
	end
end
hamming(a::AbstractStrain, b::AbstractStrain, normalization=false; ambiguous=:ignore) = hamming(a.seq, b.seq, normalization, ambiguous=ambiguous)
"""
	hamming(P::Profile{A}, Q::Profile{A}, normalization=false; 
				ambiguous=:ignore) where A <: BioSymbol
"""
function hamming(P::Profile{A}, Q::Profile{A}, normalization=false; 
				ambiguous=:ignore) where A <: BioSymbol
	length(P) != length(Q) && @error "Profiles of different lengths"
	out = 0.
	for (i,(p,q)) in enumerate(zip(P,Q))
		val = 1.
		for a in BioSequences.alphabet(A)
			if ambiguous == :ignore && isambiguous(a) 
				val -= get(p,a,0.) + get(q,a,0.) - get(p,a,0.)*get(q,a,0.)
			else
				val -= get(p,a,0.)*get(q,a,0.)
			end
		end
		out += val
	end
	if normalization
		out /= length(P)
	end
	return out
end
function hamming(P::Profile, s::BioSequence, normalization = false; ambiguous = :ignore)
	length(P) != length(s) && @error "Sequence and profile of different lengths"
	out = 0.
	for (i,(a,p)) in enumerate(zip(s,P))
		out += 1.
		for b in keys(p)
			if b == a || (ambiguous==:ignore && (isambiguous(a)||isambiguous(b)))
				out -= p[b]
			end
		end
	end
	if normalization 
		out /= length(P)
	end
	return out
end
hamming(s::BioSequence, P::Profile) = hamming(P,s)
hamming(s::Strain, P::Profile) = hamming(s.seq, P)
hamming(P::Profile, s::Strain) = hamming(s.seq, P)


"""
	consensus(P::Profile{A}) where A<:BioSymbol
	consensus(X::Array{<:BioSequence,1})

Return the consensus *sequence*, of type `LongSequence{alphabet}` where `eltype(alphabet)==A`. Does **NOT** handle ties. In case of a tie, the chosen symbol is undetermined. 

	consensus(X::Array{Strain} [, data = Dict(:strain=>"consensus")])

Return consensus of `X`, with field `data`. Does **NOT** handle ties. In case of a tie, the chosen symbol is undetermined. 
"""
function consensus(P::Profile{A}) where A<:BioSymbol
	seq = LongSequence{symbol_to_alphabet(A)}(length(P))
	for (i,p) in enumerate(P)
		if !isempty(p)
			seq[i] = findmax(p)[2]
		else
			seq[i] = ambiguous(A)
		end
	end
	return seq
end
consensus(X::Array{<:BioSequence,1}) = consensus(Profile(X))
function consensus(X::Array{<:Strain,1}, data = Dict(:strain=>"consensus"))
	p = Profile(X)
	seq = consensus(p)
	return Strain(seq, data)
end

"""
	numseq(s::LongSequence{AminoAcidAlphabet})
	numseq(S::Array{LongSequence{AminoAcidAlphabet}})
"""
function numseq(s::LongSequence{AminoAcidAlphabet})
	σ = Array{Int64,1}(undef, length(s))
	for (i,a) in enumerate(s)
		σ[i] = haskey(aa2int, a) ? aa2int[a] : 1
	end
	return σ
end
function numseq(S::Array{LongSequence{AminoAcidAlphabet}}) 
	if isempty(S)
		out = Array{Int64,2}(undef, 0, 0)
	else
		out = Array{Int64,2}(undef, length(S), length(S[1]))
		for (i,s) in enumerate(S)
			out[i,:] .= numseq(s)
		end
	end
	return out
end



