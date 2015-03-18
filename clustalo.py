import sys

from Bio.Align.Applications import ClustalOmegaCommandline

"""
The purpose of this script is to align a set of sequences using Clustal Omega.
"""

if __name__ == '__main__':
	handle = sys.argv[1]

	infile = '{0}.fasta'.format(handle)
	outfile = '{0}_aligned.fasta'.format(handle)

	cline = ClustalOmegaCommandline(infile=infile, outfile=outfile, verbose=True, auto=True, force=True)
	cline()

	