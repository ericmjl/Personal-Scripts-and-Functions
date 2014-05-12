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

def reverse_translate(seq):
    """
    This function takes in an amino acid sequence (str object) and returns the
    reverse translation of that sequence. 

    The nucleotide sequence will contain degenerate positions. 
    """
    nucleotide = ""
    for amino_acid in seq:
        nucleotide += str(dna[amino_acid])
    return nucleotide