# FluPredictibility
This repository contains the code used in the article  
``
Limited predictability of clade frequencies of seasonal influenza viruses  
Pierre Barrat-Charlaix, John Huddleston, Trevor Bedford & Richard Neher
``  

## Required libraries  
This code is entirely written in the Julia language. It uses several external modules, which can be installed by entering the following line of code in the Julia REPL (after typing `]`):   
`add BioSequences CurveFit Dates Dierckx IJulia KernelDensity LaTeXStrings Measures Plots PyPlot Random Statistics StatsBase StatsPlots`  
and  
`add https://github.com/PierreBarrat/FastaIO.jl`    
    

## Organisation  
- Folders `BioTools`, `EarthMoversDistance` and `TreeTools` contain custom Julia modules that are used to process the data and generate plots. They need to be added to variable `LOAD_PATH`, for instance by running the command `push!(LOAD_PATH, "path/to/FluPredictibility")`. They will then be usable by calling `using BioTools`.    
- Folder `data` contains trees and lists of positions in HA and NA proteins necessary for the analysis. It does *not* contain the sequence alignments used in the paper. Accession numbers of strain used are available as a supplementary file with the publication.   
- Folder `notebooks` contains jupyter notebooks that generate the figures of the paper. These notebook expect a folder `alignments` to exist, containing amino acid or nucleotide alignments to be used.   

