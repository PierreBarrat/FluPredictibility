# FluPredictibility
This repository contains the code used in the article  
``
Limited predictability of clade frequencies of seasonal influenza viruses  
Pierre Barrat-Charlaix, John Huddleston, Trevor Bedford & Richard Neher
``  

## Installing
You need Julia version 1.1 or higher to run this code, and IJulia to run the notebooks.  
If you cloned this repository to your computer, navigate to it, start a julia REPL session and type  
`]activate .`  
followed by   
`]registry add https://github.com/BioJulia/BioJuliaRegistry.git`  
and  
`]instantiate`  
This will install all dependencies. 

We will then set up IJulia with  
`]build IJulia`  
You should now be able to start jupyter notebooks with   
`using IJulia; IJulia.notebook()`  
or by launching jupyter directly from your shell. 
    

## Content
- Folder `data` contains trees and lists of positions in HA and NA proteins necessary for the analysis. It does *not* contain the sequence alignments used in the paper. Accession numbers of strain used are available as a supplementary file with the publication.   
- Folder `notebooks` contains jupyter notebooks that generate the figures of the paper. These notebook expect a folder `alignments` to exist, containing amino acid or nucleotide alignments to be used. 

