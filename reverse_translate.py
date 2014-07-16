degenerate = {
	'A' : 'GCN',
	'R' : 'MGN',
	'N' : 'AAY',
	'D' : 'GAY',
	'C' : 'TGY',
	'Q' : 'CAR',
	'E' : 'GAR',
	'G' : 'GGN',
	'H' : 'CAY',
	'I' : 'ATH',
	'L' : 'YTN',
	'K' : 'AAR',
	'M' : 'ATG',
	'F' : 'TTY',
	'P' : 'CCN',
	'S' : 'WSN',
	'T' : 'ACN',
	'W' : 'TGG',
	'Y' : 'TAY',
	'V' : 'GTN',
	'B' : 'RAY',
	'Z' : 'SAR',
	'*' : 'TAA'
	}

ecoli = {
	'A' : 'GCG',
	'R' : 'GCT',
	'N' : 'AAC',
	'D' : 'GAT',
	'C' : 'TGC',
	'Q' : 'CAG',
	'E' : 'GAA',
	'G' : 'GGC',
	'H' : 'CAT',
	'I' : 'ATC',
	'L' : 'CTG',
	'K' : 'AAA',
	'M' : 'ATG',
	'F' : 'TTT',
	'P' : 'CCG',
	'S' : 'AGC',
	'T' : 'ACC',
	'W' : 'TGG',
	'Y' : 'TAT',
	'V' : 'GTG',
	'B' : 'RAY',
	'Z' : 'SAR',
	'*' : 'TAA'
	}


codons = dict()
codons['degenerate'] = degenerate
codons['ecoli'] = ecoli


def reverse_translate(seq, organism='ecoli'):
	"""
	This function takes in an amino acid sequence as a string, and returns the 
	reverse translation of that sequence. 

	The nucleotide sequence will contain degenerate positions.
	"""
	nucleotide = ""
	for amino_acid in seq:
		nucleotide += str(codons[organism][amino_acid])
	return nucleotide