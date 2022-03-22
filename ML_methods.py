def get_short_site_data(seq, i, bwd, fwd):
  '''
  seq - str, the full DNA/RNA sequence to be sampled
    T and U will be treated as the same nucleotide
  i - int, the index of the editing site (assumed to be C)
  bwd - int, # of nts prior to editing site to be included in vector 
  fwd - int, # of nts prior to editing site to be included in vector

  Returns a list of ints, len (bwd + fwd) * 2, representing a binary vector
    Each nt represented by 2 ints in the list: is_purine, pairs_GC
      is_purine - 1 if nt is purine (AG), 0 if nt is pyrimidine (CT)
      pairs_GC - 1 if nt participates in GC pairing (GC), 0 if it doesn't (AT)
  '''
  # Create list of 0s of correct length and bound seq to be considered
  vector = [0]*((bwd + fwd + 1)*2)
  min_i = max(0, i - bwd)
  max_i = min(len(seq) - 1, i + fwd)

  sub_seq = seq[min_i:max_i + 1]

  # Changes 0 to 1 if nt belongs to correct category
  for i in range(len(sub_seq)):
    char = sub_seq[i].lower()
    if char.lower() in 'ag':
      vector[i*2] = 1
    if char.lower() in 'gc':
      vector[i*2 + 1] = 1

  return vector[:bwd*2] + vector[bwd*2 + 2:]


def gen_site_data(gene_path, fol_path):
  '''
  gene_path - str, path to the CSV list of real sites
    Must be a series of lines of format "gene,loc\n" with a header
    loc must be 1-indexed position in coding sequence
  fol_path - str, path to the folder of site files
    Filenames must be in format "gene-cds.fasta"
    Must contain files corresponding to all genes in gene_path file
  
  Returns a dict of real & fake sites in format {gene:[loc1, loc2]}
    Fake sites include 3x number of real sites, with 2 types:
      2x sites that are chosen at complete random
      1x sites that are chosen from position 25 in RNASee output
    Some genes may lack a fake site
  '''
  import rna_see, random

  # Extract real sites from gene_path file
  with open(gene_path, 'r') as f:
    lines = f.readlines()

  real_sites = {}

  for line in lines[1:]:
    gene, loc = line.split(',')
    loc = int(loc)

    if gene in real_sites:
      real_sites[gene].append(loc)
    else:
      real_sites[gene] = [loc]

  # Randomly distributes fake sites based on # of real sites
  num_real = 0
  for gene in real_sites:
    num_real += len(real_sites[gene])

  fake_sites = {}
  sites_per_gene_random = [0]*len(real_sites) # Randomly chosen sites
  sites_per_gene_see = [0]*len(real_sites) # Position 25 RNASee sites

  for i in range(num_real*2):
    sites_per_gene_random[random.randrange(len(sites_per_gene_random))] += 1
  for i in range(num_real):
    sites_per_gene_see[random.randrange(len(sites_per_gene_see))] += 1

  # Chooses random sites from list of Cs in each gene
  for gene, num in zip(real_sites, sites_per_gene_random):
    # Generates list of Cs in the file
    with open(fol_path + '/' + '%s-cds.fasta' % gene) as f:
      seq = ''.join(f.read().split('\n')[1:])

    c_pos = []
    for i in range(len(seq)):
      if seq[i].lower() == 'c':
        c_pos.append(i + 1)

    # Chooses sites
    fake_sites[gene] = []
    for _ in range(num):
      loc = random.choice(c_pos)
      while loc in real_sites[gene] or loc in fake_sites[gene]:
        loc = random.choice(c_pos)

      fake_sites[gene].append(loc)

  print('x2 Random Sites Selected\n')

  # Chooses sites similar to editing sites from position 25 RNASee
  # In case of same site being chosen already, takes a higher ranking site
  for gene, num in zip(real_sites, sites_per_gene_see):
    # Gets RNASee output of high-scoring sites
    out = rna_see.see(fol_path + '/' + '%s-cds.fasta' % gene)
    end = len(out)

    # Extracts rank-25 site
    # Chooses higher rank if site already in real or fake sites list
    for _ in range(num):
      i = 25
      loc = int(out.iloc[[i]]['pos_c'])
      while loc in real_sites[gene] or loc in fake_sites[gene]:
        i -= 1
        loc = int(out.iloc[[i]]['pos_c'])

      fake_sites[gene].append(loc)
    print('x', end='')

  print('\nx1 25+ See Sites Selected')

  return real_sites, fake_sites

def generate_vector_file(neg_file, pos_file, folder_path, out_file):
  '''
  neg_file - str, path to file of negative (non-editing) sites
    csv with header, lines in format gene,loc
  pos_file - str, path to file of positive (editing) sites
    csv with header, lines in format gene,loc
  folder_path - str, path to folder containing fasta files
    filenames should be in format gene-cds.fasta
  out_file - str, path of the output file
    
  Creates CSV file of vectors in format gene,class,vector...\n
  '''

  # Get the gene names and locations of each negative site
  with open(neg_file, 'r') as f:
    genes = f.read().split('\n')

  data = [] # Data for output to vector file

  # Generate the row data (gene name, class identifier, vector)
  for gene in genes[1:]:
    if ',' not in gene:
      continue
    
    name, loc = gene.split(',')

    with open(folder_path + '%s-cds.fasta' % name, 'r') as f:
      seq = ''.join(f.read().split('\n')[1:])

    site_info = get_short_site_data(seq, int(loc) - 1, 15, 10)
    site_info = [str(x) for x in site_info]
    data.append(name + ',0,' + ','.join(site_info) + '\n')

  # Get the gene names and locations of each positive site
  with open(pos_file, 'r') as f:
    genes = f.read().split('\n')

  # Generate the row data (gene name, class identifier, vector)
  for gene in genes[1:]:
    if ',' not in gene:
      continue
    
    name, loc = gene.split(',')

    with open(folder_path + '%s-cds.fasta' % name, 'r') as f:
      seq = ''.join(f.read().split('\n')[1:])

    site_info = get_short_site_data(seq, int(loc) - 1, 15, 10)
    site_info = [str(x) for x in site_info]
    data.append(name + ',1,'+ ','.join(site_info) + '\n')

  # Writes output file
  with open('strat_short_vectors.csv', 'w') as f:
    f.writelines(data)


from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import numpy as np
import pickle
import random

# Parameters for random_forest_creation grid search
grid_search_params = {'n_estimators':[10, 50, 100, 250, 500,1000,2000],\
                      'max_features':[3, 5, 'sqrt', 10, 20, None],\
                      'max_depth':[int(x) for x in np.linspace(10, 60, num = 6)],\
                      'min_samples_leaf':[1, 2, 4,10,20,50,100]}

def random_forest_creation(positive_vectors, negative_vectors, neg_selected):
  '''
  positive_vectors - str, path to csv with vectors of editing sites; no header
    Col 1 is gene name
    Col 2 is class identifier; 0 - non-editing, 1 - editing
    Col 3 on is the vector, eg one generated by get_short_site_data
  negative_vectors - str, path to csv file with ALL non-editing sites
    same format as positive_vectors
  neg_selected - str, path to csv file with non-editing sites for training
    format gene,loc with a header

  Creates 3 files:
    csv of vectors of test sites (with proportional pos + neg to full dataset)
    txt file of metrics of the optimized random forest model
    pickle containing optimized random forest as determined by a grid search
  '''
  # Gets list of vectors
  with open(positive_vectors, 'r') as f:
    pos_lines = f.read().strip().split('\n')

  if pos_lines[0].split(',')[1] not in ['0', '1']:
    pos_lines = pos_lines[1:]

  print('%d positive sites loaded' % len(pos_lines))

  with open(negative_vectors, 'r') as f:
    neg_lines = f.read().strip().split('\n')

  if neg_lines[0].split(',')[1] not in ['0', '1']:
    neg_lines = neg_lines[1:]

  print('%d negative sites loaded' % len(neg_lines))

  # Gets a list of non-editing sites for use in training set
  # Allows selection of a subset of non-editing sites for more even ratio
  with open(neg_selected, 'r') as f:
    strat_set = f.read().strip().split('\n')[1:]

  neg_genes = {}

  for line in strat_set:
    gene, loc = line.split(',')
    if gene in neg_genes:
      neg_genes[gene].append(int(loc))
    else:
      neg_genes[gene] = [int(loc)]

  # Creates sets of included and excluded editing and non-editing sites
  vector_set = []
  bool_set = []
  site_set = []

  for line in pos_lines: # All editing sites are included
    vector_set.append([int(x) for x in line.split(',')[2:]])
    bool_set.append(int(line.split(',')[1]))
    site_set.append(line.split(',')[0])

  nvector_set = []
  nbool_set = []
  nsite_set = []

  for line in neg_lines: # Only non-editing sites in neg_lines included
    gene, loc = line.split(',')[0].split('_')
    if gene in neg_genes and int(loc) in neg_genes[gene]:
      vector_set.append([int(x) for x in line.split(',')[2:]])
      bool_set.append(int(line.split(',')[1]))
      site_set.append(line.split(',')[0])
    else:
      nvector_set.append([int(x) for x in line.split(',')[2:]])
      nbool_set.append(int(line.split(',')[1]))
      nsite_set.append(line.split(',')[0])

  vector_array = np.array(vector_set)
  bool_array = np.array(bool_set)
  site_array = np.array(site_set)

  print('Data processed for ML; 3 arrays of following lengths')
  print(len(vector_array), len(bool_array), len(site_array))
  print()

  # Creates train and test set using a 7:3 split
  # Saves testing sites, including:
  #   All non-editing sites from neg_selected not included in training (30%)
  #   30% of non-editing sites excluded from neg_selected
  x_train, x_test, y_train, y_test, train_site, test_site = train_test_split(vector_set,\
                                                                             bool_set,\
                                                      site_set,\
                                                      test_size=0.3,\
                                                      stratify=bool_set,\
                                                      random_state=11)


  test_site = list(test_site) + random.sample(nsite_set, int(len(nsite_set)*0.3))

  # Edit test site filename here
  with open('strat_test_sites.txt', 'w') as f:
    f.write('\n'.join(test_site))

  print('%d test sites recorded' % len(test_site))

  # Creates & fits random forest classifier; prints some basic stats
  cv = StratifiedKFold(n_splits=5, shuffle=True)

  print('Running Grid Search\n')
  grid_search = GridSearchCV(estimator=RandomForestClassifier(),\
                             param_grid=grid_search_params, n_jobs=8,\
                             scoring='average_precision', cv=cv)

  grid_search.fit(x_train, y_train)
  cv_score = grid_search.best_score_
  test_score = grid_search.score(x_test, y_test)
  y_pred = grid_search.predict(x_test)

  print('CV score: %.3f\nTest score: %.3f' % (cv_score, test_score))

  # Outputs some basic metrics into a file
  scores = {}

  scores['f1'] = metrics.f1_score(y_test, y_pred)
  scores['f1_micro'] = metrics.f1_score(y_test, y_pred, average='micro')
  scores['f1_macro'] = metrics.f1_score(y_test, y_pred, average='macro')
  scores['auc_roc'] = metrics.roc_auc_score(y_test, y_pred)
  scores['recall'] = metrics.recall_score(y_test, y_pred)
  scores['precision'] = metrics.precision_score(y_test, y_pred)
  scores['confusion_matrix'] = metrics.confusion_matrix(y_test, y_pred)

  # Edit metric output filename here
  with open('strat_grid_search_results.txt', 'w') as f:
    f.write('%s\nCV score: %.3f\nTest score: %.3f\n' %\
            (str(grid_search.best_params_), cv_score, test_score))
    for key in scores:
      f.write('%s,%s\n' % (key, str(scores[key])))

  # Pickles best model; edit pickle filename here
  with open('strat_best_model.pickle', 'wb') as f:
    pickle.dump(grid_search.best_estimator_, f)

def scan_gene(ML_pickle, file_path, subset=[]):
  '''
  ML_pickle - str, filepath to a pickle file containing a ML model
    ML model must take as input the output of get_short_site_data function
  file_path - str, filepath to the CDS file for the gene in question
  subset - iterable of ints, 1-indexed pos of sites to check

  Returns list of 1-indexed pos of all predicted editing sites in given gene
  '''
  # Prepares ML model and sequence string for use
  with open(ML_pickle, 'rb') as pickle_file:
    model = pickle.load(pickle_file)

  with open(file_path, 'r') as gene_file:
    seq = ''.join(gene_file.read().split('\n')[1:])

  chosen_sites = []

  # Uses ML model to predict class of all C sites
  if subset:
    for loc in subset:
      if seq[loc].lower() != 'c':
        continue

      test_vector = [get_short_site_data(seq, loc-1, 15, 10)]

      if model.predict(test_vector)[0] and model.predict(test_vector)[0] != -1:
        chosen_sites.append(loc)
  else:
    for i in range(len(seq)):
      if seq[i].lower() != 'c':
        continue

      test_vector = [get_short_site_data(seq, i, 15, 10)]

      if model.predict(test_vector)[0] and model.predict(test_vector)[0] != -1:
        chosen_sites.append(i + 1)

  return chosen_sites
