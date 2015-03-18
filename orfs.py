"""
Author: Eric J. Ma
Date: 2015-03-17

Purpose: 
The purpose of this python module is to provide a function to find the 
longest ORF from a coding sequence.
"""

def longest_orfs(list_of_seqrecords):
    aa_sequences = []
    for i, record in enumerate(list_of_seqrecords):
        longest_protein = SeqRecord(id=record.id, seq='')
        for frame in range(3):
            length = 3 * ((len(record) - frame) // 3)
            translation = record.seq[frame:frame + length].translate()
            for pro in translation.split("*"):
                if len(pro) > len(longest_protein.seq):
                    longest_protein.seq = pro
                    
        longest_protein.seq = trim_to_start_codon(longest_protein.seq)
        aa_sequences.append(longest_protein)
        
    return aa_sequences

def trim_to_start_codon(translated_seq_obj):
    idx = str(translated_seq_obj).index('M')
    return translated_seq_obj[idx:]