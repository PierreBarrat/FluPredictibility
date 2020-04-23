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
`]build IJulia`  

Alternatively, you can type  
`]add https://github.com/PierreBarrat/FluPredictibility`  
in a julia REPL session to clone and install dependencies at the same time. Rune the IJulia relative commands after.  

You should now be able to start jupyter notebooks with   
`using IJulia; IJulia.notebook()`  
or by launching jupyter directly from your shell. 
    

## Organisation  
- Folder `data` contains trees and lists of positions in HA and NA proteins necessary for the analysis. It does *not* contain the sequence alignments used in the paper. Accession numbers of strain used are available as a supplementary file with the publication.   
- Folder `notebooks` contains jupyter notebooks that generate the figures of the paper. These notebook expect a folder `alignments` to exist, containing amino acid or nucleotide alignments to be used. 

