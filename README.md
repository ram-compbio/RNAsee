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

Scikit-learn version 0.24.2 was used to create the random forest model; while no errors have been observed when unpickling in a more recent version of scikit-learn, it is recommended that version 0.24.2 is installed. This can be done through a Python package manager like Anaconda or via pip using the command below:

```pip install scikit_learn==0.24.2```

RNAsee can be directly downloaded from this Github repository and run as invidiual scripts. On average, installation of scikit-learn and RNAsee should take less than a few minutes.

## Usage
Before using RNAsee, run the test file test.py. This will test the four major models (rules-based, random forest, union, and intersection) on predicting potential APOBEC3-mediated editing sites in 3 test files. The output should look as follows:

```
>>> import test
Running rules-based tests...
[1603 2097 1779 1462]
[136 590 438 79]
[43 2197 3006 2572 2590 1212]

Running ML tests...
[170, 325, 384, 459, 1431, 1462, 2002, 2100, 2143]
[132, 136, 649, 738]
[43, 333, 344, 1054, 1876, 2634, 2646, 2831, 3006, 3033]

Running union tests...
[170, 325, 384, 459, 1431, 1462, 2002, 2100, 2143, 1603, 2097, 1779]
[132, 136, 649, 738, 590, 438, 79]
[43, 333, 344, 1054, 1876, 2634, 2646, 2831, 3006, 3033, 2197, 2572, 2590, 1212]

Running intersection tests...
[1462]
[136]
[43, 3006]
Done
```

This test should take less than 2 minutes to run on an average computer.

It is recommended that most models are run from an interactive Python session. To run either of the non-consensus models, import the corresponding module (rules_based_methods or ML_methods) and run the scan_gene function with the filename of the correct RNA sequence. Note that ML_methods.scan_gene also requires the path to the model to be used, which is 'current_model.pickle' for the most up to date model. To run either or both consensus models, import consensus and run the union, intersection, or both_consensus functions with the filename of the correct RNA sequence. 

Coding sequence files are included in the examples and CCDS.tar.gz folders. Additional coding sequence files can be downloaded from https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi. Non-coding sequence FASTA files may also be used, but note that RNAsee's performance has not been benchmarked on full RNA sequences or DNA sequences.