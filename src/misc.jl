using DelimitedFiles
using FluPredictibility.BioTools

"""
	`true` if position only has one amino acid 
"""
function ispoly(p)
	fref = first(p.data)[2]
	fmax, aref = findmax(fref.freq)
	if fmax > 0.95
		for f in values(p.data)
			if get(f.freq, aref, 0) <= 0.95
				# println(p.i)
				# println(aref)
				# println(f.freq[aref])
				return true
			end
		end
	else
		return true
	end
	return false
end

"""
"""
function export_strains(aln::String, outfile, seqtype=:aa)
	headers = ["strain", "?", "EPI", "date", "?", "?", "?", "?","?","lab","database"]
	fp = Flu.FluPop(aln, seqtype, headers, flulineage="h3n2", segment="ha", ignore_read_errors=true)

	out = Array{String,2}(undef, length(fp.strains), 4)
	for (i,s) in enumerate(values(fp.strains))
		out[i,:] .= [s[:strain], s[:EPI], s[:lab], s[:database]]
	end

	header=["Strain" "EPI" "Lab."  "Database"];
	writedlm(outfile, vcat(header, out), "\t");
	return out
end