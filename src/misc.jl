using DelimitedFiles
using DataFrames
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

"""
	compute_activity_table(intrajectories, dt, max_date)

Return a `DataFrame` containing activity statistics of input trajectories.   
- `dt`: width of timebins trajectories are based on.
- `max_date`: Maximum date. Length of the table will be `div(max_date, dt)+1`. 
"""
function compute_activity_table(intrajectories, dt, max_date)
    imax = div(max_date, dt) + 1
    activity_table = DataFrame(:days=>collect(0:dt:max_date), :ntot=>zeros(Int64, imax),
                    :nact=>zeros(Int64, imax), :nlost=>zeros(Int64, imax), :nfixed=>zeros(Int64, imax),
                    :fact=>zeros(Float64, imax), :flost=>zeros(Float64, imax), :ffixed=>zeros(Float64, imax))

    for traj in intrajectories
        for (i,t) in enumerate([x.value for x in traj.t[1:end-1]])
            activity_table[i, :nact] += 1
            activity_table[i, :ntot] += 1
        end
        for i in (length(traj.t)):imax
            if traj.fixation == :fixed
                activity_table[i,:nfixed] += 1
                activity_table[i, :ntot] += 1
            elseif traj.fixation == :lost
                activity_table[i,:nlost] += 1
                activity_table[i, :ntot] += 1
            end
        end
    end
    activity_table[!, :fact] .= activity_table[!, :nact] ./ activity_table[!, :ntot]
    activity_table[!, :flost] .= activity_table[!, :nlost] ./ activity_table[!, :ntot]
    activity_table[!, :ffixed] .= activity_table[!, :nfixed] ./ activity_table[!, :ntot]
    return activity_table
end