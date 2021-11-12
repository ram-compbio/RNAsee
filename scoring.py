import pandas as pd
from convert import *

def check_bulge(rna,temp_i,temp_j):
    '''
    rna - str, full RNA sequence being considered
    temp_i - int, left side of stem being considered
    temp_j - int, right side of stem being considered
        should not be complementary to rna[temp_i]

    Returns bool, True if 1 nt bulge exists at temp_j
        bulge exists if rna[temp_i] is complement of rna[temp_j] + 1
    '''

    temp_j+=1
    
    if rna[temp_i] == 'c' and rna[temp_j] == 'g':
        return True
    elif rna[temp_i] == 'g' and rna[temp_j] == 'c':
        return True
    elif rna[temp_i] == 'a' and rna[temp_j] == 'u':
        return True
    elif rna[temp_i] == 'u' and rna[temp_j] == 'a':
        return True

    return False

# Calculates score of site & creates DataFrame
def palindrome(rna,loop,i,j,columns):
    '''
    rna - str, full RNA sequence being considered
    loop - int, length of loop being considered
    i - int, index of the first nt of the loop in str rna
    j - int, index of the edited C (last nt of loop) in str rna
    columns - tuple/list, names of the columns in output DataFrame
        should be 9 columns

    Calculates the score of a site given following:
        1 point per nt in the stem + 2*(number of GC pairs in stem)
        +2 points for an A or G at -1 (following the editing site)
        +2 points for a U in the loop
        -2 points for a G in the loop

    Returns empty DataFrame if site has no stem
    Returns 9-value DataFrame if site has a stem:
        temp_rna - str, the stem-loop RNA sequence
        temp_dna - str, the stem-loop sequence translated to DNA
        loop - int, the loop length
        stem - int, the stem length
        bulge - bool, whether there is a mismatch at -2 & match after
        pos_c - int, the 1-indexed position of the editing site
        stem_start - int, the index of the start of the stem in str rna
        stem_end - int, the index of the end of the stem in str rna
        strength - int, the score of the site (higher = more likely site)
    '''
    temp_i = i-1
    temp_j = j+1
    pos_c = j+1
    length = len(rna)
    stem = 0
    
    # Adds points for non-stem-based features
    strength = 2 if rna[temp_j] in 'ag' else 0
    strength += 2 if 'u' in rna[temp_i + 1:temp_j] else 0
    strength -= 2 if 'g' in rna[temp_i + 1:temp_j] else 0

    # Adds point for stem length and number of GC pairs
    bulge = False
    while temp_i >= 0 and temp_j < length:

        # 1 point for an nt pair + 2 points for G/C pair
        if rna[temp_i] == 'c' and rna[temp_j] == 'g':
            stem += 1
            strength += 3
        elif rna[temp_i] == 'g' and rna[temp_j] == 'c':
            stem += 1
            strength += 3

        # 1 point for an nt pair
        elif rna[temp_i] == 'a' and rna[temp_j] == 'u':
            stem += 1
            strength += 1
        elif rna[temp_i] == 'u' and rna[temp_j] == 'a':
            stem += 1
            strength += 1
            
        # Check for bulge at -2 position ONLY
        # Bulge exists if -2 does not match corresponding stem location,
        #   but -3 does
        # Adds correct point value
        elif temp_j - pos_c == 1 and check_bulge(rna,temp_i,temp_j):
            bulge = True
            strength += 1 if rna[temp_i] in 'aut' else 3
            temp_i -= 1
            temp_j += 2
            stem += 1
            
            continue # Continue because counters already iterated

        else:
            break # Stem is over if there is a non-bulge mismatch

        temp_i -= 1
        temp_j += 1

    if stem < 1:
        return pd.DataFrame()

    # Return to the last matching nucleotide pair
    temp_i += 1
    temp_j -= 1
    
    # Store stem-loop sequence in RNA and DNA format
    temp_rna = rna[temp_i:temp_j+1]
    temp_dna = rna2dna(temp_rna)

    # Create returned stem-loop
    df = pd.DataFrame([[temp_rna,temp_dna,loop,stem,bulge,pos_c,temp_i+1,temp_j+1,strength]], columns=columns)
    return df
