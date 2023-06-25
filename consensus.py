import rules_based_methods
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
  see_output = rules_based_methods.scan_gene(file_path, threshold=9)['pos_c']

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
  see_output = rules_based_methods.scan_gene(file_path, threshold=9)['pos_c']
  
  inter_sites = ML_methods.scan_gene(ML_model, file_path,\
                                   subset=see_output.values)

  return inter_sites

def both_consensus(file_path, ML_model='current_model.pickle'):
  '''
  file_path - str, path to FASTA file for consideration
  ML_model - str, path to pickle file with machine learning model
    Defaults to current_model.pickle in same folder

  Returns two lists of ints, union and intersection; faster than indv running
    union - list of sites considered editing by ML OR RNASee
    intersection - list of sites considered editing by ML AND RNASee
  '''
  ML_output = ML_methods.scan_gene(ML_model, file_path)
  see_output = rules_based_methods.scan_gene(file_path, threshold=9)['pos_c']

  inter_sites = []
  union_sites = []

  for loc in ML_output:
    union_sites.append(loc)
    if loc in see_output.values:
      inter_sites.append(loc)

  for loc in see_output.values:
    if loc not in union_sites:
      union_sites.append(loc)

  return ML_output, list(see_output.values), union_sites, inter_sites

class Scores:
  def __init__(self, gene, site, rules_s, ML_s):
    self.gene = gene
    self.site = site
    self.rules_s = rules_s
    self.ML_s = ML_s

  def __repr__(self):
    return '|%s_%d, %s, %s|' % (self.gene, self.site,\
                                    str(self.rules_s),\
                                    str(self.ML_s))

def both_consensus_sc(file_path, ML_model='current_model.pickle', subset=[]):
  '''
  file_path - str, path to FASTA file for consideration
  ML_model - str, path to pickle file with machine learning model
    Defaults to current_model.pickle in same folder

  Returns two lists of ints, union and intersection; faster than indv running
    union - list of sites considered editing by ML OR RNASee
    intersection - list of sites considered editing by ML AND RNASee
  '''
  ML_output = ML_methods.proba_scan_gene(ML_model, file_path, subset=subset)
  see_output = rules_based_methods.scan_gene(file_path, threshold=None,\
                                             subset=subset, make_output=False)

  gene = file_path.strip().split('\\')[-1].split('/')[-1].split('-')[0]
  gene = gene.lower()

  inter_sites = []
  union_sites = []
  union_values = []

  for loc in ML_output:
    union_values.append(loc)
    if loc in see_output['pos_c'].values:
      score = int(see_output[see_output['pos_c']==loc]['score'].values)
      scores_obj = Scores(gene, loc, score, ML_output[loc])
      inter_sites.append(scores_obj)
      union_sites.append(scores_obj)
    else:
      union_sites.append(Scores(gene, loc, None, ML_output[loc]))

  for loc in see_output['pos_c'].values:
    if loc not in union_values:
      score = int(see_output[see_output['pos_c']==loc]['score'].values)
      union_sites.append(Scores(gene, loc, score, None))

  return union_sites, inter_sites
  
