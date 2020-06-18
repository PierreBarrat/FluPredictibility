import FastaIO.writefasta
export writefasta, readfastastrains

# Closure for assigning labels to strains
let strain_number = 0
	global new_strain_number() = strain_number+=1
	global reset_strain_number() = strain_number = 0
end
function new_strain_label()
	return "strain_$(new_strain_number())"
end

"""
	writefasta(s::AbstractStrain, fields; fillvals = false)
	writefasta(f::IO, s::AbstractStrain, fields; fillvals=false)
	writefasta(S::Array{<:AbstractStrain}, fields; fillvals = false)
	writefasta(f::String, S::Array{<:AbstractStrain}, fields; fillvals=false, mode="w")

Write strain `s` to a fasta format. Header is built using `fields` from `s.data`. 
"""
function writefasta(s::AbstractStrain, fields; fillvals = false)
	# Make the header
	header = mapreduce(*, fields, init = "") do x
		if haskey(s.data, String(x))
			"$(get(s.data, String(x), 0))|"
		elseif haskey(s.data, Symbol(x))
			"$(get(s.data, Symbol(x), 0))|"
		elseif fillvals
			"?|"
		else
			@error "Field `$x` not in strain"
		end
		# "$(get(s.data, x) do ; fillvals ? "?|" : @error "Field $x not in strain"; end)|"  -- One liner, but does not handle the symbol/string 
	end[1:end-1] # Removing trailing '|'

	return (header, String(s.seq))
end
writefasta(f::IO, s::AbstractStrain, fields; fillvals=false) = writefasta(f, [writefasta(s, fields, fillvals=fillvals)])
function writefasta(f::String, S::Array{<:AbstractStrain}, fields; fillvals=false, mode="w")
	open(f, mode) do io
		for s in S
			writefasta(io, s, fields, fillvals=fillvals)
		end
	end
end
function writefasta(S::Array{<:AbstractStrain}, fields; fillvals = false)
	out = Array{Tuple{String, String},1}(undef, length(S))
	for (i,s) in enumerate(S)
		out[i] = writefasta(s, fields, fillvals = fillvals)
	end
	return out
end


function readfastastrains(f::Union{AbstractString,IO}, ::Val{A}, headerfields; separator = '|', strainfilters=[x->true], ignore_read_errors=false) where A <: BioSequences.Alphabet
	strains = Array{Strain{A},1}(undef, 0)
	nfiltered = 0
	nunread = 0
	ntot = 0
	typeof(f) <: AbstractString ? println("Reading $f...") : println("Reading alignment...")
	for (n,s) in FastaReader(f)
		dat = parse_header(n, headerfields, separator)
		st = try
			Strain(LongSequence{A}(s), dat)
		catch err
			if ignore_read_errors
				st = Strain(:aa)
			elseif occursin("X",s)
				st = Strain(LongSequence{A}(replace(s, "X"=>"N")), dat)
			elseif occursin("x",s)
				st = Strain(LongSequence{A}(replace(s, "x"=>"N")), dat)
			else
				println(s)
				error(err)
			end
		end
		if BioTools.isempty(st)
			nunread += 1
		elseif mapreduce(f->f(st), *, strainfilters, init=true)
			push!(strains, st)
		else
			nfiltered += 1
		end
		ntot += 1
	end
	println("Read $(length(strains)) strains out of $ntot. Filtered $nfiltered. Could not read $nunread")
	return strains	
end
function readfastastrains(f::Union{AbstractString,IO}, ::Val{A}, headerfields; separator = '|', strainfilters=[x->true], ignore_read_errors=false) where A <: Integer
	strains = Array{ArtificialStrain{A},1}(undef, 0)
	nfiltered = 0
	nunread = 0
	ntot = 0
	typeof(f) <: AbstractString ? println("Reading $f...") : println("Reading alignment...")
	for (n,s) in FastaReader(f)
		dat = parse_header(n, headerfields, separator)
		st = ArtificialStrain([parse(A,x) for x in s], dat)
		if BioTools.isempty(st)
			nunread += 1
		elseif mapreduce(f->f(st), *, strainfilters, init=true)
			push!(strains, st)
		else
			nfiltered += 1
		end
		ntot += 1
	end
	println("Read $(length(strains)) strains out of $ntot. Filtered $nfiltered. Could not read $nunread")
	return strains	
end
"""
	readfastastrains(f::Union{AbstractString,IO}, sequence_type::Symbol, headerfields; separator = '|', strainfilters=[x->true])
	readfastastrains(f::Union{AbstractString,IO}, sequence_type::DataType, headerfields; separator = '|', strainfilters=[x->true])

Possible symbols are `$(BioTools.sequenceymbols)`. Type of output will depend on the symbol used. 

Implemented filters are `hasdate` and `gapfilter`
"""
readfastastrains(f::Union{AbstractString,IO}, sequence_type::Symbol, headerfields; separator = '|', strainfilters=[x->true], ignore_read_errors=false) = readfastastrains(f, Val(BioTools.type(sequence_type)), headerfields, separator=separator, strainfilters=strainfilters, ignore_read_errors=ignore_read_errors)
readfastastrains(f::Union{AbstractString,IO}, sequence_type::DataType, headerfields; separator = '|', strainfilters=[x->true], ignore_read_errors=false) = readfastastrains(f, Val(sequence_type), headerfields, separator=separator, strainfilters=strainfilters, ignore_read_errors=ignore_read_errors)


"""
- `headerfields` is the list of fields that should be parsed. `?` are ignored. If the header `h` is longer than `headerfields`, the end of it is ignored. 
- `specialfields` are fields which require a special treatment. An example is the date (by default with `special_fields`), which calls `parse_date`. 
"""
function parse_header(h, headerfields, separator; specialfields=special_fields)
	sh = split(h, separator)
	dat = Dict()
	for (i,f) in enumerate(headerfields)
		if in(f, specialfields)
			dat[f] = parse_special_field[f](sh[i])
		elseif !in(f, ignored_header_fields)
			dat[f] = sh[i]
		end
	end
	return dat
end

function parse_date(s::AbstractString)
	date = missing
	try 
		date = Date(s)
	catch
		if !isnothing(match(r"-[0-9][0-9]*-XX$",s))
			# Day is missing, which does not matter
			date = Date(s[1:end-3])
		elseif !isnothing(match(r"[1-2][0-9][0-9][0-9]",s))
			# Month is missing, year is not --> date is missing
			date = missing
			# y = parse(Int64, match(r"[1-2][0-9][0-9][0-9]",hp["date"]).match)
		else
			date = missing
		end
	end
	return date
end

parse_float(s::AbstractString) = parse(Float64, s)