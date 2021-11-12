import sys
import pandas as pd

# Function to convert DNA sequence to RNA sequence
def dna2rna(dna):
    rna = []
    # clean DNA str
    dna = dna.lower().replace(' ','').replace('\n','')
    # create rna sequence with complements to dna sequence
    for i in dna:
        if i == 't':
            rna += ['u']
        elif i == 'a':
            rna += ['a']
        elif i == 'g':
            rna += ['g']
        elif i == 'c':
            rna += ['c']
        else:
            rna += [i]
            #print("Unknown nucleotide {}...exiting program".format(i))
            #sys.exit()
    rna = ''.join(rna)
    return rna


# Function to convert RNA sequence to DNA sequence
def rna2dna(rna):
    dna = []
    # clean DNA str
    rna = rna.lower().replace(' ','').replace('\n','')
    # create rna sequence with complements to dna sequence
    for i in rna:
        if i == 'u':
            dna += ['t']
        elif i == 'a':
            dna += ['a']
        elif i == 'g':
            dna += ['g']
        elif i == 'c':
            dna += ['c']
        elif i == 'n':
            dna += ['n']
        else:
            rna += [i]
            #print("Unknown nucleotide {}...exiting program".format(i))
            #sys.exit()
    dna = ''.join(dna)
    return dna
