#!/usr/bin/env python

from Bio.Seq import *
import pandas as pd
from convert import *
from scoring import *
import argparse
import os, math

def aa_change(pos_c,prot,rna,j):
    '''
    pos_c - int, 1-indexed location of editing site in RNA
    prot - str, protein sequence corresponding to RNA
    rna - str, full RNA coding sequence being considered
    j - int, 0-indexed location of editing site in str rna
    
    Returns tuple of (amino acid change, index of change in protein sequence)
        amino acid change - str, <original> -> <changed> or "synonymous"
        index of change - int
    '''
    pos = math.ceil(pos_c/3.0)
    pos = int(pos)

    # Edit the C>U
    temp_rna = list(rna)
    temp_rna[j] = 'u'
    temp_rna = ''.join(temp_rna)
    
    # Translate the edited sequence
    temp_prot = translate(temp_rna)
    if prot[pos-1] != temp_prot[pos-1]:
        f = "{} -> {}".format(prot[pos-1],temp_prot[pos-1])
    else:
        f = "synonymous"
    return f,pos

   
def see(seq, is_rna=False, cutoff=None):
    '''
    seq - str, location of fasta file to be considered
        must be a coding sequence, not full DNA/RNA, for aa_change
        must have a comment/non-code sequence on first line
    is_rna - bool, True if seq refers to a FASTA containing "U" instead of "T"
    cutoff - int, number of sites to be returned; all returned for cutoff=None

    Generates a tsv file containing sorted potential APOBEC3A/G editing sites
        in same directory as the fasta file
    Returns sorted DataFrame of potential APOBEC3A/G editing sites
        rna - str, sequence of RNA in stem-loop surrounding editing site
        dna - str, sequence of DNA corresponding to stem-loop
        loop_len - int, length of loop (3 or 4) w/C at end
        stem_len - int, length of stem surrouding loop
        bulge - bool, whether there is a mismatch at -2 position in stem
        pos_c - int, 1-indexed location of C in the full nt sequence
        begin - int, first nt in the stem-loop
        end - int, last nt in the stem-loop
        score - int, score of the site; primary sorting factor
        aa_chagne - str, amino acid change caused by C>U edit at the site
        aa_pos - int, position of altered amino acid in the protein seq
    '''
    out = seq.split("/")[-1].replace(".fasta","")
    
    # Read in input fasta
    with open(seq, 'r') as f:
        f.readline() # Skip comment on first line
        lines = f.readlines()
    
    if is_rna == True:
        # Input was RNA sequence
        rna = ''.join(lines)
        rna = rna.lower()
    else:
        # Input was DNA; gets corresonding RNA sequence
        dna = ''.join(lines)
        dna = dna.lower()
        rna = dna2rna(dna)
    
    prot = translate(rna) # corresponding protein sequence
    
    # Create dataframes
    columns=('rna','dna','loop_len','stem_len','bulge','pos_c','begin','end','score')
    df_all = pd.DataFrame(columns=columns)
    df3 = pd.DataFrame(columns=columns)
    df4 = pd.DataFrame(columns=columns)
    
    size = len(rna)
    i = 0
    j = 3

    # Iterate through all Cs within the sequence and check editing probability
    while j < size:
        if rna[j] == 'c':
            pos_c = j+1
            if (rna[j-1] == 'u' or rna[j-1] == 'c'):
                current = None
                max_stem = 0
                max_bulge = False

                # Check for stem-loop surrounding site with a 4-nt loop
                loop = 4
                temp_df = palindrome(rna,loop,i,j,columns)
                a,b = aa_change(pos_c,prot,rna,j)
                temp_df['aa_change'] = a
                temp_df['aa_pos'] = int(b)
                if not temp_df.empty:
                    df4 = df4.append(temp_df, ignore_index=True, sort=False)

                # Check for stem-loop surrounding site with a 3-nt loop
                shift_i = i+1
                loop = 3
                temp_df = palindrome(rna,loop,shift_i,j,columns)
                a,b = aa_change(pos_c,prot,rna,j)
                temp_df['aa_change'] = a
                temp_df['aa_pos'] = int(b)
                if not temp_df.empty:                    
                    df3 = df3.append(temp_df, ignore_index=True, sort=False)

        i += 1
        j += 1

    # Combine DataFrames and sort according to score
    df_all = df_all.append([df4, df3], ignore_index=True)
    df_all = df_all.sort_values(by=['score', 'stem_len'],\
                                ascending=[False, False])
    df_all = df_all.drop_duplicates(subset='pos_c', keep='first')

    # Generate file and return DataFrame
    if cutoff != None:
        df_all[:cutoff].to_csv('{}-top{}.tsv'.format(out, cutoff),sep='\t',float_format='%.2f')
    else:
        df_all.to_csv('{}.tsv'.format(out),sep='\t',float_format='%.2f')
    return df_all[:cutoff] if cutoff else df_all


if __name__ == '__main__':
    # Arugments
    parser = argparse.ArgumentParser(prog="RNAsee", description="Search a sequence for putative RNA editing sites by APOBEC3A/G.")
    parser.add_argument('sequence', metavar='sequence', type=str, help="FASTA file for DNA seqeunce input.")
    parser.add_argument('--cutoff','-c', type=int, help="Return only the top X high-scoring values.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s v0.2')
    parser.add_argument('--rna', action="store_true", help="The FASTA file is an RNA sequence (only A, U, G, C).")
    args=parser.parse_args()
    # Set variables
    seq = args.sequence
    is_rna = args.rna
    cutoff = args.cutoff
    df = see(seq, is_rna, cutoff=cutoff)
