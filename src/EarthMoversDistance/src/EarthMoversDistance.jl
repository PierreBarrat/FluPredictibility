module EarthMoversDistance

import JuMP
import COSMO, MathOptInterface
# import BioTools
using FluPredictibility.BioTools
using BioSequences: BioSequence

export EMD

"""
	EMD(f,P,Q)

Compute the earth mover's distance between distributions `P` and `Q`, using `f` as a metric. 
Enumerating `P` (resp. `Q`) should yield pairs `(x,p)` (resp. `(y,q)`) representing state `x` and its weight `p`. Distance between elements is `f(x,y)`.
"""
function EMD(f::Function,P,Q)
	# Building the distance matrix
	d = Array{Float64, 2}(undef, length(P), length(Q))
	for (i,(x,p)) in enumerate(P)
		for (j,(y,q)) in enumerate(Q)
			d[i,j] = f(x,y)
		end
	end
	# Building the JuMP model
	model = JuMP.Model(JuMP.with_optimizer(COSMO.Optimizer))
	JuMP.set_silent(model)
	JuMP.@variable(model, w[i=1:length(P), j=1:length(Q)] >= 0)
	for (i,(x,p)) in enumerate(P)
		JuMP.@constraint(model, sum(w[i,j] for j in 1:length(Q)) ==  p)
	end
	for (j,(y,q)) in enumerate(Q)
		JuMP.@constraint(model, sum(w[i,j] for i in 1:length(P)) ==  q)
	end
	JuMP.@objective(model, Min, sum(sum(w.*d)))
	JuMP.optimize!(model)
	# Looking at the solution
	if JuMP.termination_status(model) != MathOptInterface.OPTIMAL
		@warn "Termination status of optimization: $(JuMP.termination_status(model))"
	end
	flow = JuMP.value.(w)
	emd = JuMP.objective_value(model)

	return (emd=emd, flow=flow)
end

##
## Specialized functions for BioSequences and BioTools
##
"""
	EMD(f::Function, X::Array{<:BioSequence,1},Y::Array{<:BioSequence,1})	
	EMD(X::Array{<:BioSequence,1},Y::Array{<:BioSequence,1})
"""
function EMD(f::Function, X::Array{<:BioSequence,1},Y::Array{<:BioSequence,1})
	P = sequence_distribution(X)
	Q = sequence_distribution(Y)
	ret = EMD(f, P, Q)
	flow = Dict()
	for (i,(x,p)) in enumerate(P)
		for(j,(y,q)) in enumerate(Q)
			flow[x,y] = ret.flow[i,j]
		end
	end
	return (emd=ret.emd, flow=flow), P, Q
end
EMD(X::Array{<:BioSequence,1},Y::Array{<:BioSequence,1}) = EMD(BioTools.hamming, X, Y)
EMD(f::Function, X::Array{<:BioTools.AbstractStrain,1}, Y::Array{<:BioTools.AbstractStrain,1}) = EMD(f, 
	[x.seq for x in X], [y.seq for y in Y])
EMD(X::Array{<:BioTools.AbstractStrain,1}, Y::Array{<:BioTools.AbstractStrain,1}) = EMD([x.seq for x in X], [y.seq for y in Y])
EMD(P::Dict{<:BioSequence,Float64}, Q::Dict{<:BioSequence,Float64}) = EMD(BioTools.hamming, P, Q)

"""
	sequence_distribution(X::Array{<:BioSequence})
	sequence_distribution(X::Array{<:AbstractStrain})

Convenience function. Transform an array `X` of sequences or strains into a distribution form usable by the EMD functions, *i.e.* a `Dict{BioSequence => weight::Float64}`. 
"""
function sequence_distribution(X::Array{<:BioSequence})
	d = Dict{BioSequence, Float64}()
	for s in X
		d[s] = get(d, s, 0.) + 1. /length(X)
	end
	return d
end
function sequence_distribution(X::Array{<:BioTools.Strain})
	d = Dict{BioSequence, Float64}()
	for s in X
		d[s.seq] = get(d, s.seq, 0.) + 1. /length(X)
	end
	return d
end

end