import pandas as pd
from convert import *
from secondary_struct import *

# Bulge check function
def check_bulge(rna,temp_i,temp_j):
    temp_j+=1
    if rna[temp_i] == 'c' and rna[temp_j] == 'g':
        return True
    elif rna[temp_i] == 'g' and rna[temp_j] == 'c':
        return True
    elif rna[temp_i] == 'a' and rna[temp_j] == 'u':
        return True
    elif rna[temp_i] == 'u' and rna[temp_j] == 'a':
        return True
    else:
        return False

# Function to determine the length of the pallindromic sequence
# Stem length
def pallindrome(rna,best,loop,length,pos_c,i,j,columns):
    df = pd.DataFrame()
    temp_i = i-1
    temp_j = j+1
    stem = 0
    bulge = False
    while temp_i >= 0 and temp_j < length:
        if rna[temp_i] == 'c' and rna[temp_j] == 'g':
            stem += 1
        elif rna[temp_i] == 'g' and rna[temp_j] == 'c':
            stem += 1
        elif rna[temp_i] == 'a' and rna[temp_j] == 'u':
            stem += 1
        elif rna[temp_i] == 'u' and rna[temp_j] == 'a':
            stem += 1
        # Check for bulge at -2 position ONLY
        elif temp_j - pos_c == 1 and check_bulge(rna,temp_i,temp_j):
            bulge = True
            temp_i -= 1
            temp_j += 2
            stem += 1
            continue
        else:
            break
        # Keep searching...
        temp_i -= 1
        temp_j += 1
    if stem < 1:
        return df
    temp_i += 1
    temp_j -= 1
    # Temp store stem-loop sequence
    temp_rna = rna[temp_i:temp_j+1]
    # Get dna complement of rna stem-loop sequence
    temp_dna = rna2dna(temp_rna)
    # Get secondary structure and energy
    ss, mfe = get_ss(temp_rna)

    # Leave commented for now...need to test more
    '''
    # Get secondary structure and energy
    # of a larger subsequence (length of 40nt)
    if bulge:
        diff = temp_j - temp_i + 1
    else:
        diff = temp_j - temp_i
    change = (40 - diff) / 2
    if temp_i - change >= 0 and temp_j + change <= len (rna):
        temp_rna_large = rna[temp_i-change:temp_j+change+1]
        ss_large, mfe = get_ss(temp_rna_large)
    elif temp_i - change < 0:
        temp_rna_large = rna[0:temp_j+1+temp_i]
        ss_large, mfe = get_ss(temp_rna_large)
    elif temp_j + change > len(rna):
        diff = len(rna) - temp_j
        temp_rna_large = rna[temp_i-diff:len(rna)+1]
        ss_large, mfe = get_ss(temp_rna_large)
    '''

    # Store all stem-loop info
    df = pd.DataFrame([[temp_rna,temp_dna,ss,best,loop,stem,bulge,pos_c,temp_i+1,temp_j+1,mfe]], columns=columns)
    return df
