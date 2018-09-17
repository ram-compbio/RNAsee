import sys

# Function to convert DNA sequence to RNA sequence
def dna2rna(dna):
    rna = []
    # clean DNA str
    dna = dna.lower().replace(' ','').replace('\n','')
    #print(dna)
    #print(dna[::-1])
    # create rna sequence with complements to dna sequence
    for i in dna:
        if i == 't':
            rna += ['a']
        elif i == 'a':
            rna += ['u']
        elif i == 'g':
            rna += ['c']
        elif i == 'c':
            rna += ['g']
        else:
            print("Unknown nucleotide {}...exiting program".format(i))
            sys.exit()
    rna = rna[::-1]
    rna = ''.join(rna)
    return rna

def pallindrome(rna,length,i,j):
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
        temp_i -= 1
        temp_j += 1
    if stem < 4:
        return
    temp_i += 1
    temp_j -= 1
    print ''.join(rna[temp_i:temp_j+1])
    print("Stem length = {}".format(stem))
    print("End positions of stem-loop sequence are {}, {}\n".format(temp_i,temp_j))



#f = open("sdhb-partial_cds.fasta",'r')
#f = open("sdhb-complete_cds.fasta",'r')
f = open("chromosome-1.fasta",'r')
f.readline()
l = f.readlines()
dna = ''.join(l)
dna = dna.lower().replace(' ','').replace('\n','')

# Input was DNA
# Convert to RNA
rna = dna2rna(dna)
print rna

size = len(rna)
i = 0
j = i+3
while j < size:
    if rna[j] == 'c':
        #print rna[i:j+1]
#        if rna[i:j+1] == 'uauc':
        if rna[i:j+1] == 'cauc' \
                or rna[i:j+1] == 'cacc' \
                or rna[i:j+1] == 'ccuc' \
                or rna[i:j+1] == 'cuuc' \
                or rna[i:j+1] == 'uauc':
            pallindrome(rna,size,i,j)
#        elif rna[j-1] == 'u':
#            pallindrome(rna,size,i,j)
    i += 1
    j += 1
