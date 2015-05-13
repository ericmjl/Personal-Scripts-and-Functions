'''
locate_mutant: will take two sequences and locate the single position at which they differ

build_mutant_dict: this will take a list of pairs of sequences and build a dictionary with positions as 
				keys, values being dictionarys with mutation pairs i.e. (A,T) each with associated number
				of occurrences.

generate_random_graphs: this function will generate a specified number of random genotype networks for a 
				sequence of a given length

'''

import networkx as nx 
from fuzzywuzzy import fuzz
from Levenshtein import *
from Levenshtein import distance
import re
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import string
import random

def locate_mutation(seq1, seq2):
	'''
	Sequences will be nodes with an edge (meaning 1 AA mutation difference) connecting them.
	This function will locate the mutation betweent the two nodes and return the position that
	this mutation would have been in the reference sequence as well as the AA change as a tuple
	i.e. (Iso --> Leu) = (I, L)
	'''
	if distance(seq1, seq2) == 1:
		for i in range(0, max(len(seq1), len(seq2))):
			if seq1[i] != seq2[i]:
				return i, seq1[i], seq2[i]
	else:
		pass
		

	# if distance(seq1, seq2) != 1:
	# 	return seq1, seq2, " These sequences have more than one mutation!"
	# else:
	# 	for i in range(0, max(len(seq1), len(seq2))):
	# 		if seq1[i] != seq2[i]:
	# 			return i, seq1[i-10:i+10], seq2[i-10:i+10]

def determine_conserved_sites(aligned_seqs):
	sequence_length = 0
	for seq in aligned_seqs:
		if len(seq.seq) > sequence_length:
			sequence_length = len(seq.seq)
	# print sequence_length
	site_dict = {i: [] for i in range(sequence_length)}
	for seq_num in range(len(aligned_seqs)):
		for k in range(sequence_length):
			if aligned_seqs[seq_num][k] not in site_dict[k]:
				site_dict[k].append(aligned_seqs[seq_num][k])
	conserved_sites = list()
	# print site_dict
	for pos in site_dict.keys():
		# print len(site_dict[key])
		if len(site_dict[pos]) == 1:
			conserved_sites.append(pos)
	return conserved_sites

def get_blosum(aa):
	aas = 'ARNDCQEGHILKMFPSTWYVBZX*'
	blosum_dict = ''

def generative_model(aligned_seqs):
	model_size = len(seqs)
	seed_sequence = seqs[random.randint(0, model_size)]
	sequence_length = len(seed_sequence)
	model_seqs = list()
	model_seqs.append(seed_sequence)
	aas = 'ARNDCQEGHILKMFPSTWYVBZX*'
	for i in range(0, model_size):
		print 'q'






