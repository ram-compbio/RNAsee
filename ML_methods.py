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


import numpy
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
def random_forest_creation(vectors, negative_loc, positive_loc):
  '''
  vectors - str, path to a csv file with vectors for dataset; no header
    Col 1 is gene name
    Col 2 is class identifier; 0 - fake site, 1 - real site
    Col 3 on is the vector, eg one generated by get_short_site_data
  negative_loc - str, path to csv file with fake sites
    all lines should be of format gene,loc with a header row
  positive_loc - str, path to csv file with real sites
    all lines should be of format gene,loc with a header row

  Returns tuple RandomForestClassifier, list of training sites used
  '''
  # Gets list of vectors
  with open(vectors, 'r') as f:
    vector_lines = f.read().split('\n')

  if ',' not in vector_lines[-1]:
    vector_lines = vector_lines[:-1]
  if vector_lines[0].split(',')[1] not in ['0', '1']:
    vector_lines = vector_lines[1:]

  # Gets a list of negative and positive sites in same order as vector file
  # (if vector file was generated using method in this file)
  gene_loc = []
  with open(negative_loc, 'r') as f:
    for line in f.readlines()[1:]:
      if ',' in line:
        gene_loc.append(tuple(line.strip().split(',') + [0]))
        
  with open(positive_loc, 'r') as f:
    for line in f.readlines()[1:]:
      if ',' in line:
        gene_loc.append(tuple(line.strip().split(',') + [1]))

  # Separates vector lines into the vectors and class identifiers (bools)
  vector_set = []
  bool_set = []

  for line in vector_lines:
    vector_set.append([int(x) for x in line.split(',')[2:]])
    bool_set.append(int(line.split(',')[1]))

  vector_array = numpy.array(vector_set)
  bool_array = numpy.array(bool_set)

  # Creates train and test set using a 7:3 split
  # Also saves genes and locations for training set to exclude from testing
  x_train, x_test, y_train, y_test, \
           loc_train, loc_test = train_test_split(vector_array, \
                                                      bool_array,\
                                                      gene_loc,\
                                                      test_size=0.3)

  # Creates & fits random forest classifier; prints some basic stats
  clf = RandomForestClassifier()
  clf.fit(x_train, y_train)
  y_pred = clf.predict(x_test)

  print('Recall:', metrics.recall_score(y_test, y_pred))
  print('Precision:', metrics.precision_score(y_test, y_pred))
  print('Accuracy:', metrics.accuracy_score(y_test, y_pred))
  print('F1 Score:', metrics.f1_score(y_test, y_pred))

  return clf, loc_train


#alpha = svm_creation('vector_file.csv')
import pickle

def scan_gene(ML_pickle, file_path):
  '''
  ML_pickle - str, filepath to a pickle file containing a ML model
    ML model must take as input the output of get_short_site_data function
  file_path - str, filepath to the CDS file for the gene in question

  Returns list of 1-indexed pos of all predicted editing sites in given gene
  '''
  # Prepares ML model and sequence string for use
  with open(ML_pickle, 'rb') as pickle_file:
    model = pickle.load(pickle_file)

  with open(file_path, 'r') as gene_file:
    seq = ''.join(gene_file.read().split('\n')[1:])

  chosen_sites = []

  # Uses ML model to predict class of all C sites
  for i in range(len(seq)):
    if seq[i].lower() != 'c':
      continue

    test_vector = [get_short_site_data(seq, i, 15, 10)]

    if model.predict(test_vector)[0] and model.predict(test_vector)[0] != -1:
      chosen_sites.append(i + 1)

  return chosen_sites
