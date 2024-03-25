"""
	translate(infasta::String, outfasta::String; CDS = [(1,:end)])

Translate DNA sequences in fasta file `infasta` to amino acids, and write result in fasta file `outfasta`.   

Use `CDS` to indicate coding regions, *e.g.* `CDS = [(start1,stop1), (start2, stop2)]`. Default `CDS = [(1,:end)]` will use all the DNA sequence. 
"""
function translate(infasta::String, outfasta::String; CDS = [(1,:end)])
	open(outfasta, "w") do of
		for (header, seq) in FastaReader(infasta)
			s = aa""
			for (start, stop) in CDS
				if stop == :end
					s *= BioSequences.translate(LongDNASeq(replace(seq[start:end], r"[EFIJLOPQUXZ-]"i=>"N")))
				else
					s *= BioSequences.translate(LongDNASeq(replace(seq[start:stop], r"[EFIJLOPQUXZ-]"i=>"N")))
				end
			end
			writefasta(outfasta, [(header, s)], "a")
		end
	end
end

"""
	replace_conflicting_symbols(s::AbstractString, seqtype::Symbol)

Replace characters in `s` that do not match the sequence type with the ambiguous symbol. Return a new string.
"""
function replace_conflicting_symbols(s::AbstractString, seqtype::Symbol)
	if seqtype == :dna
		R = Regex("[^$(_DNAAlphabet)]")
	elseif seqtype == :aa
		R = Regex("[^$(_AAAlphabet)]")
	elseif seqtype == :rna
		R = Regex("[^$(_RNAAlphabet)]")
	else
		@error "Possible symbols: `:dna`, `:aa`, `:rna`"
	end
	return replace(s, R=>Char(ambiguous(type(seqtype))))
end

"""
	replace_conflicting_symbols(infasta::String, outfasta::String, seqtype::Symbol)
"""
function replace_conflicting_symbols(infasta::String, outfasta::String, seqtype::Symbol)
	open(outfasta, "w") do of
		for (header, seq) in FastaReader(infasta)
			s = replace_conflicting_symbols(seq, seqtype)
			writefasta(outfasta, [(header, s)], "a")
		end
	end
end
