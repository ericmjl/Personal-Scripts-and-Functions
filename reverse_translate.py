dna = {
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
    '*' : 'TRR'
    }

# reverse_translate takes in only 1 parameter: "seq", which is an amino acid sequence.
# Ensure that there is only a sequence of letters present.
def rt(seq):
    nucleotide = ""
    for amino_acid in seq:
        nucleotide += str(dna[amino_acid])
    return nucleotide