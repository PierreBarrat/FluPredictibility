## For parsing JSON into dlmreadable files
epitope_positions = []; 
epitope_perauthor = Dict(); 
for f in readdir("./")
       if !isnothing(match(r".json", f))
              a = JSON.Parser.parsefile(f)
              epitope_perauthor[a["name"]] = []
              for (n,mp) in a["map"]
                     i0 = div(FluTools.gene_positions["h3n2","ha"][n][1], 3) + 1
                     for i in keys(mp)
                            push!(epitope_positions, i0+parse(Int, i))
                            push!(epitope_perauthor[a["name"]], i0+parse(Int,i))
                     end
                     println("$(a["name"]) / $(n) : $(length(mp)) positions")
              end
       end
end
       unique!(epitope_positions)
       sort!(epitope_positions)
       for (n,v) in epitope_perauthor
       sort!(v)
       end

## Loop for trajectories
param = get_aa_parameters()
authors = ["all", "koel", "luksza", "shih", "wolf"]
files = ["h3n2/data/epitopes/ha/aligned_h3n2_ha_aa_$(a)epitopes.fasta" for a in authors]
at_epitopes = Dict()
outplot_epitopes = Dict()
for (a,f) in zip(authors, files)
       println(a)
       at_epitopes[a] = get_trajectories(f, param)
       outplot_epitopes[a] = Pfix_vs_frequency(at_epitopes[a], param)
end

