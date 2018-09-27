import pandas as pd
from convert import *
from secondary_struct import *

# Function to determine the length of the pallindromic sequence
# Stem length
def pallindrome(rna,best,loop,length,i,j,columns):
    temp_i = i-1
    temp_j = j+1
    stem = 0
    while temp_i >= 0 and temp_j < length:
        if rna[temp_i] == 'c' and rna[temp_j] == 'g':
            stem += 1
        elif rna[temp_i] == 'g' and rna[temp_j] == 'c':
            stem += 1
        elif rna[temp_i] == 'a' and rna[temp_j] == 'u':
            stem += 1
        elif rna[temp_i] == 'u' and rna[temp_j] == 'a':
            stem += 1
        else:
            break
        # Keep searching...
        temp_i -= 1
        temp_j += 1
#    if stem <= 1:
#        return
    temp_i += 1
    temp_j -= 1
    # Temp store stem-loop sequence
    temp_rna = rna[temp_i:temp_j+1]
    # Get dna complement of rna stem-loop sequence
    temp_dna = rna2dna(temp_rna)
    # Get secondary structure and energy
    ss, mfe = get_ss(temp_rna)
    # Store all stem-loop info
    df = pd.DataFrame([[temp_rna,temp_dna,ss,best,loop,stem,temp_i,temp_j,mfe]], columns=columns)
    return df
