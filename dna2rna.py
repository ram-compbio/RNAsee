import sys
import pandas as pd

# Function to convert DNA sequence to RNA sequence
def dna2rna(dna):
    rna = []
    # clean DNA str
    dna = dna.lower().replace(' ','').replace('\n','')
    print(dna)
    #print(dna[::-1])
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
            print("Unknown nucleotide {}...exiting program".format(i))
            sys.exit()
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
        else:
            print("Unknown nucleotide {}...exiting program".format(i))
            sys.exit()
    dna = ''.join(dna)
    return dna

# Function to determine the length of the pallindromic sequence
# Stem length
def pallindrome(rna,loop,length,i,j,columns):
    temp_i = i-1
    temp_j = j+1
    stem = 0
    while temp_i >= 0 and temp_j < length:
        if rna[temp_i] == 'c' and rna[temp_j] == 'g':
            stem += 1
            #print rna[temp_i:temp_j+1]
        elif rna[temp_i] == 'g' and rna[temp_j] == 'c':
            stem += 1
            #print rna[temp_i:temp_j+1]
        elif rna[temp_i] == 'a' and rna[temp_j] == 'u':
            stem += 1
            #print rna[temp_i:temp_j+1]
        elif rna[temp_i] == 'u' and rna[temp_j] == 'a':
            stem += 1
            #print rna[temp_i:temp_j+1]
        else:
            break
        # Keep searching...
        temp_i -= 1
        temp_j += 1
    if stem < 3:
        return
    temp_i += 1
    temp_j -= 1
    # Temp store stem-loop sequence
    temp_rna = rna[temp_i:temp_j+1]
    # Get dna complement of rna stem-loop sequence
    temp_dna = rna2dna(temp_rna)
    # Store all stem-loop info
    if loop == 4:
        tetra = 1
    else:
        tetra = 0
    df = pd.DataFrame([[temp_rna,temp_dna,tetra,loop,stem,temp_i,temp_j]], columns=columns)
    return df



#f = open("sdhb-wt.fasta",'r')
#f = open("sdhb-partial_cds.fasta",'r')
#f = open("pcgf3-mrna.fasta",'r')
f = open("sdhb-complete_cds.fasta",'r')
#f = open("chromosome-1.fasta",'r')
f.readline()
l = f.readlines()
dna = ''.join(l)
#dna = dna.lower().replace(' ','').replace('\n','')

# Input was DNA
# Convert to RNA
rna = dna2rna(dna)
print rna

# Create dataframe
columns=('rna','dna','tetra','loop_len','stem_len','begin','end')
df = pd.DataFrame(columns=columns)

size = len(rna)
i = 0
j = i+3
while j < size:
    if rna[j] == 'c':
        if rna[i:j+1] == 'cauc' \
                or rna[i:j+1] == 'cacc' \
                or rna[i:j+1] == 'ccuc' \
                or rna[i:j+1] == 'cuuc' \
                or rna[i:j+1] == 'uauc':
            loop = 4
            df = df.append(pallindrome(rna,loop,size,i,j,columns), ignore_index=True)
        else:
            # I should check for tri- and penta-loops here
            shift_i = i+1
            loop = 3
            df = df.append(pallindrome(rna,loop,size,shift_i,j,columns), ignore_index=True)
            shift_i = i-1
            loop = 5
            df = df.append(pallindrome(rna,loop,size,shift_i,j,columns), ignore_index=True)
    i += 1
    j += 1
# Sort df by stem length (descending)
df = df.sort_values(by=['stem_len'], ascending=False)
df = df.sort_values(by=['tetra'], ascending=False)
print df
