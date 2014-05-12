from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from collections import Counter

import os
import pandas as pd


def count_letters(seqrecord):
    counter = Counter()
    for letter in seqrecord.seq:
        counter[letter] += 1
        
    return counter

def count_fraction_N(counts):
    counter = Counter()
    fraction_N = float(counter['N']) / sum(v for k, v in counts.iteritems())
    return fraction_N

def average_quality_score(seqrecord):
    """
    This returns the average quality score across all the sequenced
    positions.
    """
    quality_scores = seqrecord.letter_annotations["phred_quality"]
    average_quality_scores = sum(i for i in quality_scores) / float(len(quality_scores))
    return average_quality_scores

# Get current working directory.
cwd = os.getcwd()

# Get all of the AB1 files present in the directory.
ab1_files = [f for f in os.listdir(cwd) if f.split('.')[1] == 'ab1']

# Data processing below.
overall_stats = []
for f in ab1_files:
    print "currently processing ab1 file: %s" % f
    handle = open(f, 'rB')
    
    sequence = SeqIO.read(handle, 'abi')
    counts = count_letters(sequence)
    fraction_n = count_fraction_N(counts)
    qs = average_quality_score(sequence)
    overall_stats.append({'sequencing result':f, 
                          'letters count':counts,
                          'fraction N':fraction_n,
                          'quality score':qs})
    print "counts: %s, fraction N: %s, quality score: %s" % (counts, fraction_n, qs)
    
    if average_quality_score(sequence) > 30:
    
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence.seq) 

        blast_record = NCBIXML.read(result_handle) 

        results = []
        for alignment in blast_record.alignments:
            results.append({'title':alignment.title, 'score':alignment.hsps[0].score})

        pd.DataFrame(results).to_csv(f + '.csv')
        
    else:
        pass

pd.DataFrame(overall_stats).to_csv('Sequencing Statistics.csv')