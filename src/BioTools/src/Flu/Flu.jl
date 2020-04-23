module Flu

using BioTools
using TreeTools
using BioSequences
using Dates, DelimitedFiles, StatsBase, JSON

import StatsBase.entropy
import Base.length

export FluPop, AAFluPop
export PosEvo
export bin_by_date!
# Should I implement a flu strain that has segment / flustrain fields and a ntseq field? 
# For now, I'll do this for AAs only. 
# One way might also be to do codon alphabets in BioSequence 
mutable struct FluPop{T<:AbstractStrain}
	strains::Dict{String, T}
	datebin::Dict{Tuple{Date,Date},Array{T,1}}
end
function FluPop(; strains = Dict{String,Strain}(), 
				datebin = Dict{Tuple{Date,Date},Array{eltype(values(strains)),1}}()) 
	return FluPop(strains, datebin)
end
function FluPop(S::Array{<:AbstractStrain,1})
	return FluPop(strains = Dict(String(x[:strain])=>x for x in S))
end
"""
	add_strain_field!(fp::FluPop, field, default=0.)
"""
function add_strain_field!(fp::FluPop, field, default=0.)
	for s in values(fp.strains)
		s.data[field] = default
	end
	nothing
end
"""
	consensus(fp::FluPop) 
"""
consensus(fp::FluPop) = consensus([x for x in values(fp.strains)])

mutable struct PosEvo{A}
	i::Int64
	alphabet::Array{A,1}
	data::Dict{Tuple{Date,Date}, SiteFrequency{A}}
end
function PosEvo(A::DataType, i::Int64) 
	return PosEvo(i, Array{A,1}(undef,0), Dict{Tuple{Date,Date}, SiteFrequency{A}}())
end


"""
	mutable struct FrequencyTraj{A}
"""
mutable struct FrequencyTraj{A}
	# All fields except `.data ` are initialized when calling the `all_trajectories` function. 
	# Info on what trajectory represents
	i::Int64
	val::A # Symbol that trajectory represents
	# Data
	date::Date # Date at t=0
	t::Array{Day,1} # Days
	freq::Array{Float64,1}
	pop::Array{Int64,1} # Size of (sampled) population at that time.
	strains::Array{<:Array{<:AbstractString, 1}, 1}
	# Fixation 
	index::Dict{Symbol, Union{Int64, Missing}} # Storing indices of start / end / activity (ie t=0)
	fixation::Symbol # :fixed, :lost, :poly
	# Misc. data
	data::Dict # Expected fields -- 
end
length(t::FrequencyTraj) = length(t.t)
# getproperty(t, field) = hasfield(t, field) ? getfield(t, field) : t.data[field]

include("IO.jl")
include("global.jl")
include("filtering.jl")
include("frequencies.jl")
include("trajectories.jl")
include("lbi.jl")
include("fitness.jl")

end