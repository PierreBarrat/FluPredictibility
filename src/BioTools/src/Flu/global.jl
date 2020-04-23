 # Flu specific
const lineages = ["h3n2" , "h1n1pdm", "vic", "yam"]
# From genebank file of strain A/Beijing/32/1992
const CDS = Dict(("h3n2","na") => [(4,1410)],
	("h3n2","ha") => [(1,48),(49,1035),(1036,1698)])
const gene_positions = Dict(("h3n2","na") => Dict("NA" => (4,1410)),
	("h3n2","ha") => Dict("SigPep" => (1,48), "HA1" => (49,1035), "HA2" => (1036,1698)))
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

outliers = Dict{Union{String,Missing},Array{String}}(missing=>String[])
for l in lineages
	p = dirname(pathof(BioTools)) * "/Flu/config/outliers_$(l).txt"
	if isfile(p)
		outliers[l] = vec(readdlm(p, String))
	else
		outliers[l] = String[]
	end
	p = dirname(pathof(BioTools)) * "/Flu/config/my_outliers_$(l).txt"
	if isfile(p)
		append!(outliers[l], vec(readdlm(p,String)))
	end
end
