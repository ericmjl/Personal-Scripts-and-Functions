"""
Author: Justin Zabilansky
Date: 9 June 2014

Purpose of this set of custom functions:

To find the best match for a small query string in a larger reference 
string and locate its position

"""


from fuzzywuzzy import fuzz


def overlapping_nmer(string, n):
    """
    This method takes in a string and breaks it into overlapping n-mer chunks as a dictionary.
    
    Input:
    - string: a string
    - n: the length of th n-mer
    
    Output:
    - string_dict: a dictionary that of the structure {'string1':position, 'string2':position+1 etc...}
    """
    string_dict = {}
    for i,letter in enumerate(string):
        if i < len(string)-n:
           string_dict[string[i:i+n]] = i
    return string_dict
    
def best_nmer(ref_str, query_str, padding=0):
    """
    best_nmer takes in a reference string and and finds the strings where it best matches
    
    Input:
    - ref_str: string to be searched
    - query_str: string that is being searched for
    - padding: make nmers being searched larger than query string if desired
    
    Output:
    - tuple with the following structure (list of best string matches, best score, dictionary of reference nmers) 
    """
    max_score = 0
    max_list = []
    ref_dictionary = overlapping_nmer(ref_str, len(query_str)+padding)
    for ref_nmer in ref_dictionary:
        score = fuzz.partial_ratio(query_str, ref_nmer)
        if score > max_score:
            max_list = []
            max_list.append(ref_nmer)
            max_score = score
        elif score == max_score:
            max_list.append(ref_nmer)
    return (max_list, max_score, ref_dictionary)
    
def find_match_position(ref_str, query_str, padding=0, shift=0):
    """
    Takes in a reference string and returns the position which matches the query string best
    
    Input:
    - ref_str: string to be searched
    - query_str: string that is being searched for
    - padding: make nmers being matched longer than the query string
    - shift: return a position within nmer instead of start of nmer
    
    Output:
    - nmer_position: list of positions of residues of interest that match the query string best 
    """
    max_list, max_score, ref_dictionary = best_nmer(ref_str, query_str, padding=padding)
    nmer_position = []
    for string in max_list:   
        nmer_position.append(ref_dictionary[string] + shift)
    return nmer_position
        
        