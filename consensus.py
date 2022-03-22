import rna_see
import ML_methods

def union(file_path, ML_model='current_model.pickle'):
  '''
  file_path - str, path to FASTA file for consideration
  ML_model - str, path to pickle file with machine learning model
    Defaults to current_model.pickle in same folder
  cutoff - int, number of top RNASee sites to consider

  Returns list of ints, sites considered editing sites by ML model OR RNASee
  '''
  union_sites = ML_methods.scan_gene(ML_model, file_path)
  see_output = rna_see.see(file_path, threshold=9)['pos_c']

  for loc in see_output.values:
    if loc not in union_sites:
      union_sites.append(loc)

  return union_sites

def intersection(file_path, ML_model='current_model.pickle'):
  '''
  file_path - str, path to FASTA file for consideration
  ML_model - str, path to pickle file with machine learning model
    Defaults to current_model.pickle in same folder

  Returns list of ints, sites considered editing sites by ML model AND RNASee
  '''
  see_output = rna_see.see(file_path, threshold=9)['pos_c']
  inter_sites = ML_methods.scan_gene(ML_model, file_path,\
                                   subset=see_output.values)

  return inter_sites

def both_consensus(file_path, ML_model='current_model.pickle'):
  '''
  file_path - str, path to FASTA file for consideration
  ML_model - str, path to pickle file with machine learning model
    Defaults to current_model.pickle in same folder
  cutoff - int, number of top RNASee sites to consider

  Returns two lists of ints, union and intersection; faster than running both
    union - list of sites considered editing by ML OR RNASee
    intersection - list of sites considered editing by ML AND RNASee
  '''
  ML_output = ML_methods.scan_gene(ML_model, file_path)
  see_output = rna_see.see(file_path, threshold=9)['pos_c']

  inter_sites = []
  union_sites = []

  for loc in ML_output:
    union_sites.append(loc)
    if loc in see_output.values:
      inter_sites.append(loc)

  for loc in see_output.values:
    if loc not in union_sites:
      union_sites.append(loc)

  return union_sites, inter_sites
  
