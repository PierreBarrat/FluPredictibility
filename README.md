# FluPredictibility
This repository contains the code used in the article  

> Limited predictability of clade frequencies of seasonal influenza viruses   
> Pierre Barrat-Charlaix, John Huddleston, Trevor Bedford & Richard Neher
  

## Installing
You need Julia version 1.10 or higher to run this code, and IJulia to run the notebooks.  
First, clone this repository to your computer. Navigate to it and start a julia session by typing `julia --project=.`. 
Alternatively, run `julia` and type `using Pkg; Pkg.activate(".")` to activate the project environment. 
Finally, run `using Pkg; Pkg.instantiate()`. This will download all the dependencies given in the `Manifest.toml` file. 

You should now be able to run the notebooks. Still in the julia session, run `using IJulia; IJulia.notebook()` (this will potentially prompt you to install jupyter through the Conda.jl package, explanations [here](https://github.com/JuliaLang/IJulia.jl)). Jupyter will start in a browser window, from which you can navigate to the FluPredictibility folder. Double-click any notebook to run it. 
    

## Content
**Important**:
- Folder `data` contains trees and lists of positions in HA and NA proteins necessary for the analysis. It does *not* contain the sequence alignments used in the paper. Accession numbers of strain used are available as a supplementary file with the publication.   
- Folder `notebooks` contains jupyter notebooks that generate the figures of the paper. These notebook expect a folder `alignments` to exist, containing amino acid or nucleotide alignments to be used. Alignments should be put there in the fasta format. 

