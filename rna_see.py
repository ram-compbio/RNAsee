import sys
import pandas as pd
from convert import *
from pallindrome import *

f = open("test/xpo1-mrna.fasta",'r')
#f = open("test/pcgf3-mrna.fasta",'r')
#f = open("test/rnh1-full.fasta",'r')
#f = open("test/sdhb-full.fasta",'r')
#f = open("test/lrp10-full.fasta",'r')
#f = open("test/app-full.fasta",'r')
#f = open("test/ap2a1-full.fasta",'r')
f.readline()
l = f.readlines()
dna = ''.join(l)

# Input was DNA
# Convert to RNA
rna = dna2rna(dna)

# Create dataframe
columns=('rna','dna','best_loop','loop_len','stem_len','begin','end')
df_all = pd.DataFrame(columns=columns)
df3 = pd.DataFrame(columns=columns)
df4 = pd.DataFrame(columns=columns)
df4_best = pd.DataFrame(columns=columns)
df5 = pd.DataFrame(columns=columns)

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
            best = 1
            df4_best = df4_best.append(pallindrome(rna,best,loop,size,i,j,columns), ignore_index=True)
        else:
            # Check other tetra-loops anyway
            best = 0
            loop = 4
            df4 = df4.append(pallindrome(rna,best,loop,size,i,j,columns), ignore_index=True)
            # Check for tri- and penta-loops here
            shift_i = i+1
            loop = 3
            df3 = df3.append(pallindrome(rna,best,loop,size,shift_i,j,columns), ignore_index=True)
            shift_i = i-1
            loop = 5
            df5 = df5.append(pallindrome(rna,best,loop,size,shift_i,j,columns), ignore_index=True)
    i += 1
    j += 1
# Sort df by stem length (descending)
df4_best = df4_best.sort_values(by=['stem_len'], ascending=False)
df4 = df4.sort_values(by=['stem_len'], ascending=False)
df3 = df3.sort_values(by=['stem_len'], ascending=False)
df5 = df5.sort_values(by=['stem_len'], ascending=False)
# Merge dfs 
df_all = df_all.append([df4_best,df4,df3,df5], ignore_index=True)
print df_all
