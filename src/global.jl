is_flu_outlier() = nothing

 # Flu specific
const lineages = ["h3n2" , "h1n1pdm", "vic", "yam"]
# From genebank file of strain A/Beijing/32/1992
const CDS = Dict(("h3n2","na") => [(4,1410)],
	("h3n2","ha") => [(1,48),(49,1035),(1036,1698)],
	("h1n1pdm", "na") => [(9,1418)],
	("h1n1pdm", "ha") => [(21,71), (72,1052), (1053,1718)])

const gene_positions = Dict(("h3n2","na") => Dict("NA" => (4,1410)),
	("h3n2","ha") => Dict("SigPep" => (1,48), "HA1" => (49,1035), "HA2" => (1036,1698)),
	("h1n1pdm", "na") => Dict("NA" => (9,1418)), 
	("h1n1pdm", "ha") => Dict("SigPep"=>(21,71), "HA1"=>(72,1052), "HA2"=>(1053, 1718)))

const substitution_rate = Dict(("h3n2","ha")=>4e-3, ("h3n2","na")=>3e-3) # Number of substitutions per site per year for NUCLEOTIDES
# Note: it's quite similar for AAs 
global coalescence_time = Year(6)
global lbi_integration_time = Day(100)
# Settings from the seasonal-flu nextstrain module
const segments = ["ha","na"]
const nextstrain_regions = ["north_america", "south_america", "europe", "china", "oceania","southeast_asia", "japan_korea", "south_asia", "west_asia", "africa"]
# 
const flu_usual_header_fields = ["strain", "virus", "", "date", "region", "country", "", "", "", "segment"]
const augur_all_header_fields = ["strain", "virus", "isolate_id", "date", "region", "country", "division", "location", "passage", "authors", "age", "gender"]
flu_usual_filters(flulineage) = [!is_flu_outlier(flulineage), BioTools.hasdate, s->BioTools.gapfilter(s,threshold=0.05)]
#
function get_standard_coordinates(i::Int64; lineage="h3n2", gene="ha")
	idx = sort([(k,x[2]/3) for (k,x) in Flu.gene_positions["h3n2", "ha"]], by=x->x[2])
	offset=0
	for (k,(ref,ix)) in enumerate(idx)
		if i <= ix
			return (ref, Int64(i-offset))
		end
		offset = ix
	end
	error("AA position $i too large for gene ($lineage,$gene)")
end

outliers = Dict{Union{String,Missing},Array{String}}(missing=>String[])
for l in lineages
	p = dirname(pathof(FluPredictibility)) * "/config/outliers_$(l).txt"
	if isfile(p)
		outliers[l] = vec(readdlm(p, String))
	else
		outliers[l] = String[]
	end
	p = dirname(pathof(FluPredictibility)) * "/config/my_outliers_$(l).txt"
	if isfile(p)
		append!(outliers[l], vec(readdlm(p,String)))
	end
end
