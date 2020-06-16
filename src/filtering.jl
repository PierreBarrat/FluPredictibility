"""
	bin_by_date!(fp::FluPop; start=:auto, last=:auto, binwidth=Day(121), binspacing=Day(121))

Bin `fp` by dates, from `start` to `last`. 
"""
function bin_by_date!(fp::FluPop{S}; 
	start=:auto,
	stop=:auto,
	binwidth=Day(121),
	binspacing=Day(121)) where S
	
	fp.datebin = Dict{Tuple{Date,Date}, Array{S,1}}()
	# Start and end dates
	if start == :auto
		startdate = findmin([x[:date] for x in values(fp.strains)])[1]
	else
		startdate = Date(start)
	end
	if stop == :auto
		stopdate = findmax([x[:date] for x in values(fp.strains)])[1]
	else
		stopdate = Date(stop)
	end
	# 
	now = startdate + binwidth
	while now < stopdate # If dlat - dstart != 0 mod[binwidth], the last bin will be smaller than others and is not considered here
		fp.datebin[now-binwidth, now] = Array{S,1}(undef, 0)
		now += binspacing
	end

	binlist = sort(collect(keys(fp.datebin)))
	for s in values(fp.strains)
		date = s[:date]
		if date < binlist[end][2] && date > binlist[1][1]
			db = _find_datebin(date, binlist)
			push!(fp.datebin[db], s)
		end
	end	

	nothing
end

"""
Assume sorted datebins
"""
function _find_datebin(date::Date, datebins::Array{Tuple{Date,Date}})
	N = length(datebins) # Max
	i0 = 1 # Min
	i = div(N,2)
	found = false
	out = datebins[i]
	c = 0
	while !found
		c += 1
		if date >= datebins[i][1]
			if date < datebins[i][2]
				found = true
				out = datebins[i]
			else
				i0 = i
				i = i0 + div(N-i0,2) + 1
			end
		else
			N = i
			i = i0 + div(N-i0,2)
		end
	end
	return out
end

"""
	filter_by_region!(fp::FluPop, r::Array{<:AbstractString,1})
	filter_by_region!(fp::FluPop, r::AbstractString)
"""
function filter_by_region!(fp::FluPop, r::Array{<:AbstractString,1})
	for S in values(fp.datebin)
		idx = findall(s-> !in(s[:region], r), S)
		deleteat!(S, idx)
	end
	for s in values(fp.strains)
		if !in(s[:region], r)
			# println(s[:region])
			# println(r[1])
			# @show s[:region] == r[1]
			delete!(fp.strains, s[:strain])
		end
	end
end
filter_by_region!(fp::FluPop, r::AbstractString) = filter_by_region!(fp, [r])
"""
	filter_by_region(fp::FluPop, r)
"""
filter_by_region(fp::FluPop, r) = begin out = deepcopy(fp); filter_by_region!(out, r); return out end

"""
	filter_by_country!(fp::FluPop, r::Array{<:AbstractString,1})
	filter_by_country(fp::FluPop, r::AbstractString)
"""
function filter_by_country!(fp::FluPop, r::Array{<:AbstractString,1})
	for S in values(fp.datebin)
		idx = findall(s-> !in(s[:country], r), S)
		for i in idx
			delete!(fp.strains, S[i][:strain])
		end
		deleteat!(S, idx)
	end
end
filter_by_country(fp::FluPop, r::AbstractString) = filter_by_country!(fp, [r])
"""
	filter_by_country(fp::FluPop, r)
"""
filter_by_country(fp::FluPop, r) = begin out = deepcopy(fp); filter_by_country!(out, r); return out end

"""
	find_strains(S::Array{<:AbstractStrain}, i::Int64, val; seqtype = :aa)

Find strains in `S` that carry value `val` at position `i`. 
"""
function find_strains(S::Array{<:AbstractStrain}, i::Int64, val; seqtype = :aa)
	return S[findall(s->s.seq[i] == val, S)]
end

"""
	datebin_to_date(d)
"""
function datebin_to_date(d)
	if length(d) != 2 || d[2] < d[1]
		@error "Invalid date bin"
	end
	return d[1] + div(d[2] - d[1], 2)
end
"""
	date_to_datebin(date, binwidth)
"""
function date_to_datebin(date, binwidth)
	return (date-binwidth, date+binwidth)
end
"""
	find_datebin(date::Date, datebins)
	find_datebin(date::Date, fp::FluPop)
"""
function find_datebin(date::Date, datebins)
	found = false
	out = first(datebins)
	for db in datebins
		if db[1] <= date && db[2] > date
			out = db
			found = true
			break
		end
	end
	if !found
		@warn "date $date was not found"
	end
	return out
end
find_datebin(date::Date, fp::FluPop) = find_datebin(date, keys(fp.datebin))


function get_regions(S::Array{<:Strain,1})
	out = Dict{String,Int64}()
	for s in S
		if !in(s[:region], ["?",'?']) 
			out[s[:region]] = get(out, s[:region], 0) + 1
		end
	end
	return out
end

function numdate(d::Date)
	return year(d) + (month(d)-1) /12 + day(d) / 365
end

# Filters for IO
"""
"""
is_flu_outlier(x::AbstractStrain, flulineage) = in(x[:strain], outliers[flulineage])
function is_flu_outlier(flulineage)
	return x->is_flu_outlier(x,flulineage)
end

"""
	has_mutation(x::AbstractStrain, mutation)
	has_mutation(mutation)

Mutation `m` is expected to be of the form `(m[1], m[2])` (position, symbol). 
"""
has_mutation(x::AbstractStrain, mutation) = (x.seq[mutation[1]] == mutation[2])
function has_mutation(mutation)
	return x->has_mutation(x, mutation)
end
