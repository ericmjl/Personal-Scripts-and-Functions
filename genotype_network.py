import networkx as nx
import sys

from Bio import AlignIO
from itertools import combinations
from Levenshtein import distance

"""
Author: Eric J. Ma
Date: 2015-03-17

This simple program generates a genotype network from a multiple sequence
alignment of isolates.
"""
if __name__ == '__main__':
	filename = sys.argv[1]
	aln_type = sys.argv[2]
	g_handle = sys.argv[3]

	G = nx.Graph()

	alignment = [s for s in AlignIO.read(filename, aln_type)]

	for s1, s2 in combinations(alignment, 2):
		s1seq = str(s1.seq)
		s2seq = str(s2.seq)

		if distance(s1seq, s2seq) == 1:
			if s1seq not in G.nodes():
				G.add_node(s1seq, accessions=[s1.id])
			if s1seq in G.nodes():
				G.node[s1seq]['accessions'].append(s1.id)
			if s2seq not in G.nodes():
				G.add_node(s2seq, accessions=[s2.id])
			if s2seq in G.nodes():
				G.node[s2seq]['accessions'].append(s2.id)

			G.add_edge(s1seq, s2seq)

	nx.write_gpickle(G, g_handle)
