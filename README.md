[![GitHub release](https://img.shields.io/github/tag-pre/zfalls/RNAsee.svg)](https://github.com/zfalls/RNAsee/releases)
[![Build Status](https://travis-ci.org/zfalls/RNAsee.png)](https://travis-ci.org/zmfalls/RNAsee)
[![codecov](https://codecov.io/gh/zfalls/RNAsee/branch/master/graph/badge.svg)](https://codecov.io/gh/zmfalls/RNAsee)
<!--[![Github All Releases](https://img.shields.io/github/downloads/zfalls/RNAsee/total.svg)](https://github.com/zfalls/RNAsee/downloads/total)-->

# RNAsee
**RNA** **s**ite **e**diting **e**valuation (**RNAsee**) is a bioinformatic python script that predicts APOBEC3A/G mediated C>U edit sites for any gene.

## System Requirements
RNAsee has been tested on Windows 10 version 10.0.19045 and Linux CentOS version 7.9.2009 machines. Most modern computers should be able to run RNAsee. Must have Python 3, Biopython, scikit-learn installed.

## Installation
Python 3, Biopython, and scikit-learn are required to run RNAsee. Any version of Biopython is acceptable, and it may be installed through pip as follows:

```pip install Biopython```

Scikit-learn version 1.0.2 was used to create the random forest model; while no errors have been observed when unpickling in a more recent version of scikit-learn, it is recommended that version 1.0.2 is installed. This can be done through a Python package manager like Anaconda or via pip using the command below:

```pip install scikit-learn==1.0.2```

RNAsee can be directly downloaded from this Github repository and run as invidiual scripts. To use the ML methods or consensus methods, the current_model.tar.gz file must be uncompressed and the original current_model.pickle file must be extracted. On average, installation of scikit-learn and RNAsee should take less than a few minutes.

## Usage
Before using RNAsee, run the test file test.py. This will test the four major models (rules-based, random forest, union, and intersection) on predicting potential APOBEC3-mediated editing sites in 3 test files. The output should look as follows:

```
>>> import test
Running rules-based tests...
[1603 2097 1779 1462]
[136 590 438 79]
[43 2197 3006 2572 2590 1212]

Running ML tests...
[170, 445, 459, 1431, 1999, 2002, 2100, 2143]
[136, 649, 738]
[43, 344, 450, 1054, 1399, 1797, 1876, 2634, 2646, 2746, 2882, 2923, 3006]

Running union tests...
[170, 445, 459, 1431, 1999, 2002, 2100, 2143, 1603, 2097, 1779, 1462]
[136, 649, 738, 590, 438, 79]
[43, 344, 450, 1054, 1399, 1797, 1876, 2634, 2646, 2746, 2882, 2923, 3006, 2197, 2572, 2590, 1212]

Running intersection tests...
[]
[136]
[43, 3006]
Done
```

This test should take about 5 minutes to run on an average computer, with the ML and union tests taking the bulk of the time.

It is recommended that most models are run from an interactive Python session. To run either of the non-consensus models, import the corresponding module (rules_based_methods or ML_methods) and run the scan_gene function with the filename of the correct RNA sequence. Note that ML_methods.scan_gene also requires the path to the model to be used, which is 'current_model.pickle' for the most up to date model. To run either or both consensus models, import consensus and run the union, intersection, or both_consensus functions with the filename of the correct RNA sequence. 

Below is an example of the commands you could run to assess the cytosines in the gene SDHB using all four models, plus the expected output. Run these commands in an interactive session in the RNAsee-master folder:

```
>>> import rules_based_methods
>>> rules_based_methods.scan_gene('test/sdhb-cds.fasta') # The output of this command will be a dataframe, and a file will also be generated.
               rna             dna loop_len  ... score   aa_change aa_pos
7   ccaucuaucgaugg  ccatctatcgatgg        4  ...    15      R -> *   46.0
50   agcugccccagcu   agctgccccagct        3  ...    11      P -> L  197.0
18   gcaacuucuaugc   gcaacttctatgc        4  ...    10  synonymous  146.0
3        ccucccgag       cctcccgag        4  ...    10      R -> *   27.0
>>> rules_based_methods.scan_gene('test/sdhb-cds.fasta')['pos_c'] # To just get the site numbers as output, use this command.
7     136
50    590
18    438
3      79
Name: pos_c, dtype: object
>>>
>>> import ML_methods
>>> ML_methods.scan_gene('current_model.pickle', 'test/sdhb-cds.fasta') # The ML methods need to have a model pickle specified in addition to the FASTA file.
[136, 649, 738]
>>> 
>>> import consensus
>>> consensus.union('test/sdhb-cds.fasta')
[136, 649, 738, 590, 438, 79]
>>> consensus.intersection('test/sdhb-cds.fasta')
[136]
```

Coding sequence files are included in the examples and CCDS.tar.gz folders. Additional coding sequence files can be downloaded from https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi. Non-coding sequence FASTA files may also be used, but note that RNAsee's performance has not been benchmarked on full RNA sequences or DNA sequences.