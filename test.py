#!/usr/bin/env python
test_files = ['test/app-cds.fasta', 'test/sdhb-cds.fasta', 'test/xpo1-cds.fasta']

import rules_based_methods

print('Running rules-based tests...')

for file in test_files:
  print(rules_based_methods.scan_gene(file, threshold=9,
                                      make_output=False)['pos_c'].values)

import ML_methods

print('\nRunning ML tests...')

for file in test_files:
  print(ML_methods.scan_gene('current_model.pickle', file))

import consensus

print('\nRunning union tests...')

for file in test_files:
  print(consensus.union(file))

print('\nRunning intersection tests...')

for file in test_files:
  print(consensus.intersection(file))

print("Done")
