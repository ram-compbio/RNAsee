import requests, os
from bs4 import BeautifulSoup

URL = 'https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=GENE&DATA=%s&ORGANISM=9606&BUILDS=CURRENTBUILDS'
CCDS_URL = 'https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=%s&ORGANISM=9606&BUILDS=CURRENTBUILDS'
LINE_LEN = 70

def download_cds(gene_name, directory=''):
  '''
  gene_name - str, name of gene to be found
  Creates cds file for gene if found on CCDS site
  '''
  req  = requests.get(URL % gene_name)
  soup = BeautifulSoup(req.text, 'html.parser')

  text = soup.get_text()
  if 'Nucleotide Sequence' not in text:
    if 'Public\nHomo sapiens' in text:
      lines = text.split('\n')
      CCDS_num = ''
      for i in range(len(lines) - 3):
        if lines[i + 2] == 'Public' and lines[i + 3] == 'Homo sapiens':
          CCDS_num = lines[i]
          break
        
    else:
      print('Gene %s not found' % gene_name)
      return

    if CCDS_num:
      req = requests.get(CCDS_URL % CCDS_num)
      soup = BeautifulSoup(req.text, 'html.parser')

      text = soup.get_text()
      
  for line in text.split('\n'):
    if 'Nucleotide Sequence' not in line:
      continue

    try:
      seq = line.split('Nucleotide Sequence')[1].split('Translation')[0]
      seq = seq.split()[2].strip()

    except Exception as e:
      print(e)
      continue

    if directory:
      file_name = '%s\\%s-cds.fasta' % (directory, gene_name)
    else:
      file_name = '%s-cds.fasta' % (gene_name)
      
    with open(file_name, 'w') as file:
      file.write('> %s cds\n' % gene_name)
      for i in range(0, len(seq) - LINE_LEN, LINE_LEN):
        file.write('%s\n' % seq[i:i+LINE_LEN])
        
      file.write(seq[i+LINE_LEN:])

    print('Wrote to %s' % file_name)
    return

  print('Gene %s not found' % gene_name)
  return

import csv

class Gene_Ref:
  def __init__(self, name, chromo, pos, base):
    '''
    name - str, name of gene, lowercase
    chromo - str, chromosome on which gene is found
    pos - int, integer position of the key base
    base - str, single character base at pos
    '''
    self.name = name
    self.chromo = chromo
    self.pos = [pos]
    self.base = [base]

  def add_pos(self, pos, base):
    self.pos.append(pos)
    self.base.append(base)

  def __repr__(self):
    return 'Gene_Ref : %s | %2s | %s | %s' % (self.name, self.chromo,\
                                     self.pos, self.base)
                                   
def generate_gene_dict(csv_path):
  '''
  csv_path - path to the csv file
  returns a dict of gene_name : Gene_Ref
    can be used to access chromosome, positions, & bases for each gene
  '''
  gene_dict = {}
  
  with open(csv_path, 'r', encoding='utf-8-sig') as csv_file:
    dict_reader = csv.DictReader(csv_file)
    for row in dict_reader:
      if row['Region'] == 'exonic':
        key = row['Gene'].lower()
        
        if key in gene_dict:
          gene_dict[key].add_pos(int(row['Position']), row['Reference base'])

        else:
          ref = Gene_Ref(row['Gene'].lower(), row['Chromosome'],\
                         int(row['Position']), row['Reference base'])
          gene_dict[key] = ref

  return gene_dict

def get_chromo_seqs(genome_path):
  seqs = {}
  files = os.listdir(genome_path)

  for file in files:
    if '.fasta' not in file:
      continue

    key = file.replace('.fasta', '')

    with open('%s\\%s' % (genome_path, file), 'r') as chr_file:
      lines = chr_file.readlines()

    lines = lines[1:]
    lines = [x.strip() for x in lines]

    seqs[key] = ''.join(lines)
    print('%s finished' % file)

  return seqs

COMP_BASES = {'c':'g', 'g':'c', 't':'a', 'a':'t'}

def get_complement(seq):
  comp_seq = ''
  for base in seq[::-1]:
    comp_seq += COMP_BASES[base.lower()]
  return comp_seq

def get_pos(gene, chromo_seq, gene_dict=None, csv_path=None):
  '''
  gene - str, name of gene
    gene-cds.fasta file should be in same directory as this is run from
  chromo_seq - str, sequence of whole chromosome
  gene_dict (optional) - dict, lowercase gene names as keys
    value should be a Gene_Ref object corresponding to that base
  csv_path (optional) - str, path to csv file
    should contain fields Chromosome, Position, Reference base, and Gene
    should also contain Region (exonic, intronic, etc)
    will not be used if gene_dict is passed
  will return tuple of all positions noted in gene_dict
    relative to cds file
  '''
  if gene_dict == None and csv_path == None:
    raise ValueError('Must pass a gene_dict or csv_path to get_pos')
    
  if gene_dict == None:
    gene_dict = generate_gene_dict(csv_path)

  with open('%s-cds.fasta' % gene, 'r') as cd_file:
    lines = cd_file.readlines()[1:]
    cd_seq = ''.join([x.strip().lower() for x in lines])
    
  ref = gene_dict[gene]
  pos_list = []

  for i in range(len(ref.pos)):
    pos = ref.pos[i] - 1
    
    bases = chromo_seq['chromo%s' % ref.chromo][pos - 10:pos+11].lower()
    
    if ref.base[i] == 'G':
      bases = get_complement(bases)

    if bases[10] != 'c':
      print('Wrong base (%s, not c)' % bases[10])

    if cd_seq.count(bases) == 0:
      late_base = bases[10:]
      early_base = bases[:11]
      if cd_seq.count(late_base) == 1:
        pos_list.append(cd_seq.find(late_base) + 1)
      elif cd_seq.count(early_base) == 1:
        pos_list.append(cd_seq.find(early_base) + 11)
      else:
        print('Gene %s pos %d not added' % (gene, ref.pos[i]))
      continue

    if cd_seq.count(bases) == 1:

      pos_list.append(cd_seq.find(bases) + 11)
    else:
      print('Gene %s pos %d not added' % (gene, ref.pos[i]))

  for pos in pos_list:
    if cd_seq[pos - 1].lower() != 'c':
      print('Wrong base at finish; gene %s and pos %d' % (gene, pos))

  return pos_list

'''
if __name__ == "__main__":      
  seqs = get_chromo_seqs('..\\Full_Genome')
  gene_dict = generate_gene_dict('..\\Table_S1_exonic.csv')

  rows = []

  try:
    for gene in gene_dict:
      download_cds(gene)
      if '%s-cds.fasta' % gene not in os.listdir():
        continue
      
      positions = get_pos(gene, seqs, gene_dict)
      for pos in positions:
        rows.append({'Gene': gene, 'Position': pos})
  except:
    raise
  finally:
    with open('out_file.csv', 'w', newline='') as out_file:
      writer = csv.DictWriter(out_file, fieldnames=['Gene', 'Position'])

      writer.writeheader()
      for row in rows:
        writer.writerow(row)

  print(len(rows))
'''
missing_genes = ['aaas', 'aars1', 'abcc6', 'abcf1', 'abhd12b', 'abhd16b', 'abhd2', 'abhd3', 'abi3bp', 'ablim1', 'abr', 'acadm', 'ace2', 'acox1', 'acox3', 'acp3', 'acp5', 'actn4', 'ada2', 'adam12', 'adam17', 'adam20', 'adam29', 'adamts1', 'adamts18', 'adamts3', 'adamts5', 'adamts9', 'adamtsl1', 'adcy3', 'adcy8', 'adcy9', 'adgra1', 'adgra2', 'adgra3', 'adgrb3', 'adgre1', 'adgre2', 'adgrg1', 'adgrg2', 'adgrg4', 'adgrg6', 'adgrg7', 'adgrl3', 'adora3', 'adprs', 'adra1b', 'adra2b', 'adrb1', 'adss1', 'adss2', 'afm', 'afp', 'ago1', 'agps', 'agrp', 'agtr1', 'ahdc1', 'ahsa2p', 'aif1', 'aifm1', 'aip', 'aire', 'ak1', 'ak7', 'akap1', 'akap12', 'akap5', 'akr1c8p', 'aldh1a1', 'aldh1a2', 'aldh1b1', 'aldh3a2', 'aldh4a1', 'aldh5a1', 'aldh6a1', 'aldh8a1', 'aldh9a1', 'aldob', 'alox12', 'alox12b', 'alox15', 'alox15b', 'alox5', 'aloxe3', 'alx4', 'amfr', 'ammecr1', 'amn', 'ampd1', 'ampd2', 'anapc5', 'angpt1', 'ankib1', 'ankle2', 'ankrd10', 'ankrd18a', 'ankrd2', 'ankrd26', 'ankrd26p1', 'ankrd6', 'anpep', 'anxa8l1', 'aoc2', 'ap1g1', 'ap1s3', 'ap2a1', 'ap2s1', 'ap3b1', 'ap3d1', 'ap4e1', 'ap5z1', 'apba3', 'apeh', 'apex2', 'apoa5', 'apob', 'apoc1', 'apof', 'aprg1', 'aqp6', 'araf', 'arcn1', 'areg', 'arfrp1', 'arg1', 'arhgap12', 'arhgap31', 'arhgap45', 'arhgef10', 'arhgef28', 'arhgef4', 'arhgef40', 'arhgef9', 'arms2', 'arnt', 'arnt2', 'arsd', 'arsf', 'arsl', 'arx', 'asb11', 'asb15', 'asb16', 'asb16-as1', 'asb2', 'ascl1', 'ascl3', 'asic2', 'asip', 'asns', 'astl', 'asxl3', 'atcay', 'atf5', 'atl1', 'atoh7', 'atp5mj', 'auh', 'auts2', 'avpr1a', 'avpr2', 'axin1', 'axl', 'azin2', 'baalc', 'bace1', 'bag5', 'baiap2', 'baiap3', 'bap1', 'barx1', 'basp1', 'baz1a', 'baz2b', 'bbs1', 'bbs7', 'bbs9', 'bcar3', 'bcat2', 'bccip', 'bche', 'bcl2', 'bcl2l10', 'bcl7a', 'bdkrb1', 'bdkrb2', 'bfsp2', 'bicdl2', 'bin2', 'birc8', 'blk', 'blm', 'bmp2', 'bmp5', 'borcs5', 'bpnt2', 'bptf', 'brawnin', 'brd2', 'brme1', 'btc', 'btg1', 'btn2a2', 'btn2a3p', 'btn3a1', 'btn3a2', 'bub1b', 'bysl', 'c11orf40', 'c15orf54', 'c1orf147', 'c1orf167', 'c1orf220', 'c2', 'c3p1', 'c4bpa', 'c4orf54', 'c6', 'c7', 'c8b', 'c8g', 'c8orf17', 'c9', 'ca6', 'ca8', 'cabp2', 'cabp4', 'cacna2d1', 'cacna2d4', 'cacnb1', 'cacnb3', 'cacng1', 'cacng2', 'cacng6', 'cad', 'calm1', 'calm3', 'camk2d', 'camkk1', 'cap2', 'capn1', 'capn11', 'capn12', 'capn3', 'capn6', 'capns1', 'caps2', 'card14', 'card16', 'card6', 'cars1', 'casc2', 'cbr1', 'cbr3', 'cbwd5', 'cby2', 'ccdc137', 'ccl23', 'ccnb3', 'ccne2', 'ccnk', 'ccnyl2', 'ccr3', 'ccr4', 'ccr7', 'ccr9', 'ccrl2', 'cct4', 'cct5', 'cct7', 'cct8', 'cct8l1p', 'cct8l2', 'cd163', 'cd19', 'cd200', 'cd27', 'cd2ap', 'cd79b', 'cd83', 'cd96', 'cdc23', 'cdc25c', 'cdc27', 'cdc5l', 'cdc6', 'cdin1', 'cdk10', 'cdk15', 'cdk2', 'cdk5', 'cdk6', 'cdk8', 'cdk9', 'cdkl1', 'cdkn3', 'cds1', 'cdsn', 'cdx2', 'ceacam18', 'ceacam6', 'ceacam8', 'cefip', 'celf2', 'cenatac', 'cenpc', 'cenpe', 'cenph', 'cep43', 'cerkl', 'cers1', 'cert1', 'cfap20dc', 'cfap221', 'cfap251', 'cfap74', 'cfap91', 'chd6', 'chfr', 'chia', 'chmp4b', 'chrdl2', 'chrna1', 'chrna10', 'chrna2', 'chrna3', 'chrna4', 'chrna5', 'chrna6', 'chrna9', 'chrnb1', 'chrnb3', 'chrnb4', 'chrnd', 'chrne', 'cib1', 'cib4', 'cibar2', 'cic', 'cilk1', 'ckap2', 'ckap4', 'ckb', 'ckm', 'clcn2', 'clcn4', 'clcn5', 'clcnka', 'clic5', 'clic6', 'clk1', 'clk3', 'clk4', 'cllu1', 'clmn', 'cln3', 'clrn1', 'clstn1', 'clybl', 'cmtm2', 'cmtm6', 'cnmd', 'cnot2', 'cnp', 'cntn4', 'cntnap2', 'cntnap5', 'cog1', 'cog3', 'cog4', 'cog5', 'cog6', 'cog7', 'col10a1', 'col11a2', 'col17a1', 'col18a1', 'col19a1', 'col21a1', 'col4a1', 'col4a5', 'col4a6', 'col6a2', 'col9a1', 'colec10', 'colec12', 'coq6', 'corin', 'cplx1', 'cpne3', 'cpne7', 'cracd', 'cracdl', 'crat', 'crhr2', 'crlf1', 'crlf3', 'crnkl1', 'crtap', 'cryaa', 'cryab', 'crybg1', 'crym', 'cryz', 'csf2ra', 'csn1s1', 'csrnp2', 'cst8', 'cstl1', 'ctage1', 'ctdsp1', 'ctrl', 'ctsl3p', 'cul3', 'custos', 'cxadr', 'cxcl11', 'cxcl3', 'cxcr3', 'cxcr6', 'cxorf1', 'cyb5a', 'cyb5r3', 'cybb', 'cycs', 'cyld', 'cyp11a1', 'cyp11b1', 'cyp11b2', 'cyp19a1', 'cyp24a1', 'cyp27a1', 'cyp2d6', 'cyp2d7', 'cyp2f1', 'cyp2s1', 'cyp39a1', 'cyp3a4', 'cyp3a43', 'cyp3a5', 'cyp3a7', 'cyp4f12', 'cyp51a1', 'cyp7b1', 'cyrib', 'cyth4', 'daam2', 'dact2', 'dag1', 'dars1', 'dcaf11', 'dctn4', 'ddc', 'ddn', 'ddo', 'ddx27', 'ddx55', 'decr1', 'degs1', 'dek', 'delec1', 'dennd2b', 'dgcr8', 'dgka', 'dgkg', 'dgkq', 'dguok', 'dhrs4', 'dhrs4l2', 'dhx8', 'diaph2', 'dio2', 'dkk2', 'dlc1', 'dlgap1', 'dlk1', 'dll3', 'dll4', 'dlst', 'dlx5', 'dmbx1', 'dmp1', 'dmrta1', 'dna2', 'dnaaf10', 'dnaaf11', 'dnaaf2', 'dnaaf6', 'dnaaf8', 'dnaaf9', 'dnah10', 'dnai3', 'dnai4', 'dnai7', 'dnaja3', 'dnajb11', 'dnajb12', 'dnajb2', 'dnajb6', 'dnajb9', 'dnajc30', 'dnajc6', 'dnase1', 'dnase1l1', 'dnase2', 'dnm1p46', 'dock8-as1', 'dpcd', 'dpf1', 'dpf3', 'dpp3', 'dpp4', 'dpy19l1', 'drd5', 'dsc1', 'dsc3', 'dscr9', 'dtx2', 'dusp13', 'dusp19', 'dusp29', 'dusp6', 'dut', 'dxo', 'dync2i1', 'dync2i2', 'dynlt2', 'dynlt5', 'dyrk2', 'dyrk4', 'dytn', 'e2f1', 'e2f3', 'ebi3', 'ebna1bp2', 'ech1', 'echs1', 'edar', 'eddm3a', 'eddm3b', 'ednra', 'ednrb', 'eef1akmt4', 'efemp1', 'efs', 'egf', 'egfem1p', 'egln3', 'ehd2', 'ehhadh', 'eif2ak4', 'eif2b4', 'eif2b5', 'eif3b', 'eif3d', 'eif5a2', 'elac1', 'elac2', 'elapor1', 'elapor2', 'eloa3dp', 'eloa3p', 'elovl2', 'elovl4', 'elp1', 'elspbp1', 'emd', 'eml5', 'eng', 'eno2', 'eno3', 'enpep', 'enpp1', 'enpp2', 'enpp3', 'enpp4', 'enpp5', 'entpd6', 'epb41l2', 'epb41l4a', 'epc1', 'epcam', 'epha1', 'epha4', 'epha7', 'ephb6', 'epm2a', 'epor', 'eprs1', 'eps8', 'eral1', 'ercc5', 'ereg', 'ern2', 'erv3-1', 'ervk-18', 'ervw-1', 'esd', 'esr2', 'esrrb', 'ess2', 'etfb', 'etfdh', 'eurl', 'evc', 'evl', 'exd2', 'exoc3-as1', 'eya1', 'f10', 'f11', 'f2', 'f8', 'f9', 'fads3', 'fah', 'fam120aos', 'fam13a', 'fam13b', 'fam13c', 'fam166c', 'fam167a', 'fam205bp', 'fam209a', 'fam209b', 'fam220bp', 'fam25a', 'fam3a', 'fam71e1', 'fam71f2', 'fam72b', 'fam74a7', 'fam83c', 'fam83g', 'fam86c1p', 'fanca', 'fancd2', 'fance', 'fap', 'farp2', 'fas', 'faslg', 'fbln5', 'fbn1', 'fbrsl1', 'fbxl17', 'fbxl3', 'fbxl4', 'fbxo10', 'fbxo11', 'fbxo24', 'fbxo30', 'fbxo32', 'fbxo34', 'fbxo4', 'fbxo42', 'fbxo5', 'fbxo7', 'fbxw7', 'fdft1', 'fdx2', 'fech', 'fem1a', 'fer1l4', 'fer1l6-as2', 'ferd3l', 'fez2', 'fgf14', 'fgf20', 'fgf21', 'fgf23', 'fgf5', 'fgf8', 'fgf9', 'fgfr1', 'fgfr2', 'fgfr4', 'fhip1a', 'fhip1b', 'fhip2a', 'fhip2b', 'fhl2', 'fkbp4', 'fkbp9p1', 'fkbpl', 'flna', 'flnc', 'flot1', 'flrt2', 'flt3', 'flt4', 'fmo6p', 'fnbp1', 'fnbp4', 'foxa1', 'foxa2', 'foxa3', 'foxc1', 'foxe1', 'foxg1', 'foxh1', 'foxi1', 'foxl2', 'foxn1', 'fpgs', 'fpgt', 'fras1', 'frmd6-as1', 'fryl', 'fshr', 'fth1p19', 'ftsj3', 'fxr1', 'fxr2', 'fyn', 'g6pc1', 'gaa', 'gabbr1', 'gabbr2', 'gabpb1', 'gabra2', 'gabra4', 'gabra5', 'gabra6', 'gabrb1', 'gabrb2', 'gabrb3', 'gabre', 'gabrg2', 'gabrg3', 'gabrr1', 'gad1', 'galc', 'galk1', 'galk2', 'galns', 'galr1', 'galr3', 'gan', 'gars1', 'gas8', 'gas8-as1', 'gask1b', 'gast', 'gatb', 'gba3', 'gbe1', 'gbf1', 'gbx1', 'gc', 'gcdh', 'gch1', 'gclc', 'gcm2', 'gcna', 'gdf3', 'gdf5', 'gemin5', 'gemin6', 'get1', 'gfpt1', 'gga2', 'ggh', 'ghrh', 'ghrhr', 'gipr', 'gjb2', 'gjb6', 'gjb7', 'gla', 'gli4', 'glp1r', 'glp2r', 'glra1', 'glrx5', 'gls', 'gmfg', 'gmpr', 'gnal', 'gnb1l', 'gnb5', 'gnpat', 'gnptab', 'gnrhr', 'gns', 'golga5', 'golim4', 'gorasp2', 'gosr2', 'got2', 'gp1ba', 'gp1bb', 'gp6', 'gp9', 'gpd1', 'gpd2', 'gpi', 'gpm6a', 'gpr101', 'gpr107', 'gpr108', 'gpr15', 'gpr152', 'gpr153', 'gpr179', 'gpr32', 'gpr39', 'gpr62', 'gpr63', 'gpr65', 'gpr68', 'gpr75', 'gpr78', 'gpr84', 'gpr85', 'gpr87', 'gprc5d', 'gprc6a', 'gpt', 'gpt2', 'gpx2', 'gpx3', 'gpx5', 'gpx6', 'gria1', 'gria2', 'gria4', 'grid2', 'grik1', 'grik2', 'grin2a', 'grin2b', 'grin2c', 'grin2d', 'grin3a', 'grin3b', 'grip1', 'grm3', 'grm4', 'grm6', 'grm8', 'grxcr2', 'gsdmc', 'gsk3a', 'gspt2', 'gsr', 'gss', 'gsta1', 'gsta3', 'gtf3c3', 'guca1a', 'gucy1a1', 'gucy1b2', 'h1-1', 'h1-2', 'h1-3', 'h1-4', 'h1-5', 'h1-6', 'h1-7', 'h1-8', 'h2ac21', 'h2ap', 'h2az2', 'h2bc13', 'h2bc14', 'h2bc4', 'h2bw1', 'h3-3a', 'h3-5', 'h3c1', 'h3c15', 'h4c1', 'haao', 'hal', 'hars1', 'haus5', 'hbs1l', 'hcar1', 'hccs', 'hck', 'hcp5', 'hcrtr2', 'hdac2', 'hdac7', 'hdac9', 'hdc', 'hebp1', 'hebp2', 'hectd1', 'hectd2', 'helb', 'herc2p3', 'herc6', 'herpud1', 'hesx1', 'hey2', 'hfe', 'hgd', 'hgfac', 'hhip', 'hibch', 'hif3a', 'hivep2', 'hla-a', 'hla-b', 'hla-c', 'hla-dma', 'hla-dmb', 'hla-doa', 'hla-dob', 'hla-dpa1', 'hla-dpb1', 'hla-dqa1', 'hla-dqa2', 'hla-dqb1', 'hla-dqb2', 'hla-dra', 'hla-drb1', 'hla-drb3', 'hla-drb4', 'hla-drb5', 'hla-e', 'hla-f', 'hlcs', 'hltf', 'hmgcr', 'hmgn2p46', 'hmox1', 'hmox2', 'hnf4a', 'hnmt', 'hnrnpf', 'hnrnph3', 'hnrnpul1', 'homer2', 'homer3', 'hoxa1', 'hoxa13', 'hoxa2', 'hoxa3', 'hoxa4', 'hoxa7', 'hoxb1', 'hoxb3', 'hoxb7', 'hoxc11', 'hoxc13', 'hoxc4', 'hoxc9', 'hoxd1', 'hoxd10', 'hoxd11', 'hoxd13', 'hoxd3', 'hoxd4', 'hpgd', 'hps1', 'hr', 'hrh1', 'hrh2', 'hrh4', 'hrob', 'hsd17b12', 'hsd3b7', 'hspa1a', 'hspa1b', 'hspa1l', 'hspa2', 'hspa5', 'hspa8', 'hspb1', 'hspb2', 'hspb3', 'hyal1', 'hydin', 'iars1', 'ibsp', 'icam1', 'icam2', 'icam3', 'idh2', 'idh3b', 'ids', 'ier3', 'ier5l', 'ifnar1', 'ifnar2', 'ifng', 'ifngr1', 'ifngr2', 'igbp1', 'igf2-as', 'igfals', 'igfbp6', 'igfbp7', 'igha1', 'igha2', 'ighe', 'ighg1', 'ighg2', 'ighg3', 'ighm', 'ighv1-2', 'ighv1-46', 'ighv3-20', 'ighv3-23', 'igkc', 'igkv1-5', 'igkv1d-16', 'iglc3', 'iglc7', 'igsf1', 'igsf11', 'igsf6', 'iho1', 'ikbkb', 'ikbke', 'il10ra', 'il10rb', 'il12rb1', 'il15ra', 'il1r1', 'il22ra1', 'il22ra2', 'il2rb', 'il2rg', 'il3ra', 'il4i1', 'il4r', 'il5ra', 'il7r', 'ilf3', 'impa2', 'impdh2', 'impg2', 'inha', 'inhba', 'inpp1', 'insm2', 'insr', 'ints1', 'iqank1', 'iqcn', 'irag1', 'irag2', 'irak1', 'irak2', 'irak3', 'irak4', 'irf9', 'irx3', 'irx5', 'ism1', 'itga2b', 'itga4', 'itga5', 'itgal', 'itgb4', 'itgb7', 'itk', 'ivd', 'jak2', 'jam2', 'jph1', 'jph2', 'jph3', 'kars1', 'kash5', 'katnb1', 'katnip', 'kcne1', 'kcne2', 'kcnmb1', 'kcnmb3', 'kctd8', 'kdelr3', 'kdr', 'kdsr', 'khk', 'khnyn', 'kiaa0087', 'kiaa1522', 'kiaa1549', 'kics2', 'kif13a', 'kif13b', 'kif15', 'kif1b', 'kif1c', 'kif22', 'kif23', 'kif24', 'kif25', 'kif25-as1', 'kif3c', 'kif4b', 'kif9', 'kifbp', 'kir2dl3', 'kir3dl2', 'kl', 'klf10', 'klf2', 'klf3', 'klf5', 'klf6', 'klhdc10', 'klhl28', 'klhl3', 'klhl38', 'klhl5', 'klhl8', 'kmt2c', 'kmt2d', 'kmt2e', 'kntc1', 'kremen1', 'krit1', 'krt20', 'krt23', 'krt5', 'krtap1-3', 'krtap10-1', 'krtap10-10', 'krtap10-11', 'krtap10-12', 'krtap10-2', 'krtap10-3', 'krtap10-4', 'krtap10-5', 'krtap10-6', 'krtap10-7', 'krtap10-8', 'krtap10-9', 'krtap11-1', 'krtap12-2', 'krtap12-3', 'krtap13-2', 'krtap13-4', 'krtap15-1', 'krtap19-2', 'krtap19-4', 'krtap19-8', 'krtap20-1', 'krtap20-2', 'krtap21-1', 'krtap22-1', 'krtap25-1', 'krtap26-1', 'krtap27-1', 'krtap3-2', 'krtap4-1', 'krtap4-2', 'krtap4-3', 'krtap4-4', 'krtap4-5', 'krtap4-7', 'krtap5-3', 'krtap5-6', 'krtap5-8', 'krtap5-9', 'krtap6-3', 'krtap9-2', 'krtap9-4', 'krtap9-6', 'krtap9-7', 'krtap9-9', 'ktn1', 'lair1', 'lair2', 'lama2', 'lama3', 'lamp1', 'lamp2', 'lamp3', 'lars1', 'las2', 'lats2', 'lcat', 'ldb3', 'ldlr', 'lef1', 'leg1', 'lemd3', 'lep', 'lgals16', 'lgi2', 'lgi3', 'lgi4', 'lhx3', 'liat1', 'lilra2', 'lilra3', 'lilrb1', 'lilrb2', 'lilrb4', 'lilrb5', 'linc00114', 'linc00271', 'linc00303', 'linc00322', 'linc00336', 'linc00482', 'linc00518', 'linc00523', 'linc00526', 'linc00574', 'linc01554', 'linc01555', 'linc01556', 'linc01559', 'linc01587', 'linc01600', 'linc02694', 'linc02875', 'linc02908', 'lipc', 'lipf', 'liph', 'llgl1', 'lman1l', 'lmod1', 'lmx1b', 'lnx2', 'loricrin', 'lox', 'loxl2', 'lpl', 'lpo', 'lrba', 'lrp10', 'lrp3', 'lrp6', 'lrpap1', 'lrrc15', 'lrrc2', 'lrrc38', 'lrrc74b', 'lrtm2', 'lss', 'lta', 'ltb', 'ltb4r', 'ltbr', 'ltk', 'ly6g5b', 'ly6g6c', 'ly6g6d', 'ly75', 'ly86', 'lypla1', 'lyst', 'mab21l4', 'madd', 'maf', 'mafb', 'mageb16', 'maml2', 'man2a1', 'man2b1', 'man2c1', 'map2k1', 'map2k2', 'map2k4', 'map2k7', 'map3k11', 'map3k4', 'map3k5', 'map3k9', 'map7d3', 'mapk11', 'mapk13', 'mapk3', 'mapk8', 'mapk9', 'marchf10', 'marchf2', 'marchf3', 'marchf6', 'marchf7', 'marchf8', 'marchf9', 'mark1', 'mars1', 'mast4', 'matn2', 'matn4', 'mbd5', 'mboat1', 'mboat4', 'mbtps2', 'mc2r', 'mc3r', 'mcam', 'mccc2', 'mccd1', 'mcm2', 'mcm3', 'mcm5', 'mcm7', 'mdga1', 'mdn1', 'me3', 'mea1', 'mecom', 'med13l', 'med23', 'mef2a', 'men1', 'meox1', 'mepe', 'mertk', 'met', 'mfsd6', 'mgam', 'mical3', 'mid1', 'mideas', 'minpp1', 'mipep', 'mir1-1hg', 'mir1-1hg-a', 'mir1915hg', 'mir7-3hg', 'mis18bp1', 'misp3', 'mkks', 'mlf2', 'mln', 'mmaa', 'mmab', 'mmp10', 'mmp2', 'mmp20', 'mmp21', 'mmp7', 'mmp9', 'mms19', 'mmut', 'mnx1', 'mocos', 'mok', 'morc1', 'mpdu1', 'mpg', 'mphosph9', 'mpi', 'mpo', 'mpp4', 'mras', 'mrc2', 'mrgprg-as1', 'mrpl18', 'mrpl2', 'mrpl22', 'mrpl27', 'mrpl28', 'mrpl30', 'mrps12', 'mrps16', 'mrps18b', 'mrps2', 'mrps31', 'mrs2', 'ms4a12', 'ms4a14', 'ms4a15', 'ms4a2', 'ms4a4a', 'ms4a5', 'ms4a6a', 'ms4a6e', 'ms4a7', 'ms4a8', 'msh2', 'msh5', 'msln', 'mslnl', 'msmb', 'msr1', 'mst1r', 'msx1', 'msx2', 'mt-atp6', 'mt-atp8', 'mt-co1', 'mt-co2', 'mt-co3', 'mt-cyb', 'mt-nd1', 'mt-nd2', 'mt-nd3', 'mt-nd4', 'mt-nd4l', 'mt-nd5', 'mt-nd6', 'mt-rnr1', 'mta1', 'mtarc1', 'mtarc2', 'mtdh', 'mthfd1', 'mtnr1a', 'mtnr1b', 'mtrfr', 'mtrr', 'muc3b', 'musk', 'mx1', 'mxd3', 'myb', 'mybpc1', 'mybpc3', 'myc', 'myeov', 'myg1', 'myh10', 'myh6', 'myl6', 'myo15b', 'myo1c', 'myo1d', 'myo1f', 'myo1g', 'myo3b', 'myo5b', 'myoz2', 'naa11', 'naaa', 'naga', 'nagk', 'nanogp8', 'nanos1', 'nanos2', 'nbeap1', 'nbpf7', 'ncoa2', 'ncr1', 'ncr2', 'ncr3', 'ndfip2', 'ndrg1', 'ndrg2', 'necab3', 'nedd4l', 'nek1', 'nek3', 'nek8', 'neu1', 'neurl2', 'neurod4', 'nfe2l3', 'nfkb2', 'nfkbie', 'nfx1', 'ngfr', 'ngly1', 'ngrn', 'nhlrc1', 'niban1', 'nid2', 'nipsnap1', 'nkapl', 'nkx2-1', 'nkx2-5', 'nkx2-6', 'nkx2-8', 'nkx3-1', 'nkx6-2', 'nln', 'nmbr', 'nmnat1', 'nmrk2', 'nmu', 'nod1', 'nol3', 'nomo1', 'nop9', 'nopchap1', 'nos1', 'nosip', 'notch4', 'nox3', 'noxa1', 'npas2', 'nphp4', 'nphs1', 'nphs2', 'npr2', 'npy1r', 'nqo1', 'nqo2', 'nr0b2', 'nr1d2', 'nr1h2', 'nr1h3', 'nr1i2', 'nr1i3', 'nr2f2', 'nr3c1', 'nr3c2', 'nr4a1', 'nrg1', 'nrl', 'nrp1', 'nrp2', 'nsun5p2', 'nt5c', 'nt5e', 'ntaq1', 'ntf3', 'ntn1', 'ntrk2', 'ntrk3', 'ntsr1', 'ntsr2', 'nudt11', 'nudt5', 'nudt7', 'numa1', 'numb', 'nup153', 'nup155', 'nup188', 'nup205', 'nup58', 'nynrin', 'oas1', 'oas3', 'odad1', 'odad2', 'odad3', 'odad4', 'odf1', 'odf2', 'odf2l', 'odf4', 'ogdh', 'ogfr', 'ogt', 'olfm2', 'onecut1', 'opa3', 'oprk1', 'or10a7', 'or10ad1', 'or10c1', 'or10g2', 'or10g3', 'or10h1', 'or10h2', 'or10h3', 'or10h4', 'or10j3', 'or10p1', 'or11a1', 'or11g2', 'or11h4', 'or11h6', 'or12d2', 'or12d3', 'or13c3', 'or13c8', 'or13c9', 'or13d1', 'or13f1', 'or13h1', 'or13j1', 'or14j1', 'or1c1', 'or1e2', 'or1f1', 'or1g1', 'or1i1', 'or1j1', 'or1j2', 'or1k1', 'or1l1', 'or1l3', 'or1l4', 'or1l6', 'or1l8', 'or1n1', 'or1n2', 'or1q1', 'or2a12', 'or2a14', 'or2a2', 'or2a25', 'or2a5', 'or2ae1', 'or2b2', 'or2b6', 'or2c1', 'or2f1', 'or2f2', 'or2g2', 'or2g3', 'or2h1', 'or2h2', 'or2j1', 'or2j2', 'or2j3', 'or2l8', 'or2t1', 'or2v2', 'or2w1', 'or2z1', 'or3a1', 'or3a2', 'or3a3', 'or4d1', 'or4d2', 'or4e1', 'or4e2', 'or4k1', 'or4k13', 'or4k14', 'or4k15', 'or4k2', 'or4l1', 'or4m1', 'or4m2', 'or4n2', 'or4n4', 'or4n5', 'or4q3', 'or51j1', 'or52a4p', 'or5ac2', 'or5au1', 'or5h14', 'or5h15', 'or5h2', 'or5h6', 'or5k3', 'or5k4', 'or5v1', 'or6b1', 'or6b3', 'or6c3', 'or6c4', 'or6c6', 'or6c65', 'or6c70', 'or6c74', 'or6c75', 'or6k3', 'or6n2', 'or6s1', 'or6v1', 'or7c1', 'or7c2', 'or7d2', 'or7d4', 'or7e24', 'or7g1', 'or7g2', 'or7g3', 'or8s1', 'or8u3', 'or9a2', 'or9k2', 'orc3', 'orc6', 'orm1', 'osbp', 'osbpl10', 'osbpl3', 'osbpl5', 'osbpl9', 'otc', 'otof', 'otog', 'otor', 'otud3', 'ovos2', 'oxtr', 'p2rx2', 'p2ry10', 'p2ry12', 'p3h4', 'pabpc1l', 'pacsin1', 'padi3', 'page1', 'pah', 'pak4', 'pak5', 'pak6', 'pak6-as1', 'pank4', 'panx3', 'pappa2', 'papss1', 'papss2', 'parp2', 'pax1', 'pax2', 'pax3', 'pax4', 'pax5', 'pax8', 'pax9', 'pbx4', 'pcca', 'pccb', 'pck1', 'pcmt1', 'pcna', 'pdcd11', 'pdcd6ip', 'pde10a', 'pde11a', 'pde1c', 'pde3a', 'pde4a', 'pde4c', 'pde4d', 'pde6a', 'pde6g', 'pde7a', 'pdha1', 'pdha2', 'pdhb', 'pdlim1', 'pdlim5', 'pdss1', 'pdxdc2p', 'pdyn', 'pdzd2', 'peak1', 'peds1', 'peg3', 'penk', 'per1', 'pex14', 'pex3', 'pfas', 'pfkl', 'pgam2', 'pgap2', 'pgap6', 'pgbd1', 'pgbd3', 'pgk1', 'phf1', 'phf3', 'phip', 'phka1', 'phkb', 'phldb3', 'phox2b', 'phyh', 'pi4kap2', 'pigl', 'pik3ca', 'pik3cb', 'pim1', 'pim2', 'pisd', 'pitrm1', 'pitx1', 'pitx2', 'pitx3', 'piwil1', 'pkd1', 'pkd1l1-as1', 'pkd1l2', 'pkm', 'pkn1', 'pknox1', 'pknox2', 'pla1a', 'pla2g1b', 'pla2g2d', 'pla2g4c', 'pla2g7', 'pla2r1', 'plaat4', 'plagl1', 'plb1', 'plcb1', 'plcb2', 'plcb3', 'plcb4', 'plcd1', 'plcg2', 'plcl2', 'plek2', 'plekha4', 'plekhg1', 'plk4', 'plod3', 'plp1', 'plpp2', 'pls3', 'plxnb3', 'plxnc1', 'pnkd', 'poglut2', 'pold1', 'polr1d', 'polr1g', 'polr1h', 'polr2j2', 'polr2m', 'pom121', 'pot1', 'poteb3', 'pou2f3', 'pou3f3', 'ppib', 'ppig', 'ppil1', 'ppil3', 'ppil6', 'ppp1cb', 'ppp1cc', 'prcd', 'prdm1', 'prdm14', 'prdm5', 'prdm7', 'prdm8', 'prep', 'prickle1', 'prkacg', 'prkag3', 'prkar1a', 'prkar2b', 'prkca', 'prkcd', 'prkce', 'prkcsh', 'prkg1', 'prkn', 'prlr', 'prm3', 'prmt9', 'prnt', 'procr', 'prokr1', 'prox2', 'prpf31', 'prpf4b', 'prpf8', 'prps1', 'prrc2b', 'prss1', 'prss12', 'prss16', 'prss50', 'psapl1', 'psg1', 'psg8', 'psmd12', 'psors1c2', 'pstpip2', 'ptdss1', 'ptger1', 'ptges2', 'ptgir', 'ptgis', 'ptgr1', 'ptgs1', 'ptk6', 'ptpdc1', 'ptpn1', 'ptpn11', 'ptpn12', 'ptpn13', 'ptpn18', 'ptpn21', 'ptpn3', 'ptpn4', 'ptpra', 'ptprb', 'ptprd', 'ptprf', 'ptprh', 'ptprj', 'ptprm', 'ptprn', 'ptprn2', 'ptprs', 'ptprt', 'ptpru', 'pttg2', 'pum1', 'pycr1', 'pyy3', 'qars1', 'qpct', 'rab17', 'rab23', 'rab27a', 'rab2b', 'rab38', 'rab40a', 'rab7a', 'rac2', 'rac3', 'rad21l1', 'rad51', 'rad51b', 'raet1g', 'raly', 'rap1gap2', 'rapsn', 'rarb', 'rarres1', 'rars1', 'rassf1', 'rassf3', 'rassf9', 'rax', 'rbm23', 'rcan3', 'rdh12', 'rdh8', 'reck', 'relb', 'reps1', 'rere', 'retnlb', 'rexo4', 'rgl2', 'rgs11', 'rhag', 'rhbdd1', 'rhobtb2', 'rif1', 'rin2', 'rin3', 'ring1', 'ripk2', 'ripor2', 'ripor3', 'rln1', 'rnase13', 'rnase2', 'rnase4', 'rnase6', 'rnaseh1', 'rnf10', 'rnf8', 'rock2', 'ror2', 'rorb', 'rorc', 'ros1', 'rp1l1', 'rpia', 'rpl36al', 'rpl6', 'rpp38-dt', 'rps6kl1', 'rreb1', 'rrm2b', 'rs1', 'rsph1', 'rsph6a', 'rtn1', 'rtn4', 'rtn4ip1', 'rtn4r', 'runx1', 'runx1t1', 'rusc1', 'rusc1-as1', 'rusf1', 'ryr1', 'ryr3', 'sall4', 'samsn1', 'sardh', 'sarg', 'sars1', 'sart1', 'sart3', 'sav1', 'sbf1', 'scand2p', 'scarf2', 'scgn', 'scn1b', 'scn4b', 'scn7a', 'scnn1a', 'scnn1b', 'sctr', 'sdr39u1', 'sec23b', 'sec23ip', 'sec24b', 'sec24c', 'sec24d', 'selenbp1', 'selplg', 'sema5a', 'serpinb1', 'serpinb6', 'serpini1', 'sertad4-as', 'sesn2', 'sez6', 'sez6l', 'sfrp2', 'sftpd', 'sgcb', 'sgcd', 'sgce', 'sgk1', 'sgk2', 'sgk3', 'sh3bp4', 'sh3glb2', 'shpk', 'sirt3', 'sirt5', 'sirt6', 'six1', 'six2', 'six4', 'six6', 'six6os1', 'skap2', 'sla', 'slamf1', 'slc10a1', 'slc10a2', 'slc11a2', 'slc12a1', 'slc12a3', 'slc12a5', 'slc13a3', 'slc14a2', 'slc16a10', 'slc16a2', 'slc17a1', 'slc17a3', 'slc17a4', 'slc17a5', 'slc17a7', 'slc18a1', 'slc19a1', 'slc1a1', 'slc1a5', 'slc20a2', 'slc22a1', 'slc22a20p', 'slc22a7', 'slc25a14', 'slc25a22', 'slc26a1', 'slc26a10', 'slc26a11', 'slc26a7', 'slc26a8', 'slc27a4', 'slc28a1', 'slc29a1', 'slc2a10', 'slc2a4', 'slc2a9', 'slc30a3', 'slc30a9', 'slc35b2', 'slc35d2', 'slc36a2', 'slc36a4', 'slc37a1', 'slc38a2', 'slc38a4', 'slc38a6', 'slc39a11', 'slc39a3', 'slc39a4', 'slc39a6', 'slc40a1', 'slc44a2', 'slc44a4', 'slc45a1', 'slc4a11', 'slc4a5', 'slc5a2', 'slc5a3', 'slc5a5', 'slc5a7', 'slc6a15', 'slc6a2', 'slc6a8', 'slc7a7', 'slc7a9', 'slc8a2', 'slc8a3', 'slc9a2', 'slc9a3r1', 'slf2', 'slit2', 'smad5-as1', 'smarcal1', 'smarcc1', 'smarcc2', 'smc2', 'smc4', 'smchd1', 'smg6', 'smim23', 'smoc1', 'smpd2', 'smrp1', 'snai1', 'snca', 'sncaip', 'sncg', 'snd1-it1', 'snrpb', 'sntb2', 'snx10', 'snx15', 'snx18', 'snx30', 'snx7', 'snx8', 'soat2', 'socs1', 'socs2', 'socs3', 'sod3', 'sorbs1', 'sord', 'sos2', 'sox17', 'sox18', 'sox3', 'sp1', 'spaca1', 'spag11b', 'spag6', 'spam1', 'spanxd', 'sparcl1', 'spart', 'spata25', 'spata7', 'spg11', 'sphk2', 'spint2', 'spn', 'spock2', 'spock3', 'spon2', 'spp1', 'spr', 'spred1', 'spring1', 'sprr2b', 'sprr2d', 'spry2', 'spry4', 'sptb', 'sptbn1', 'sptbn4', 'sqor', 'sqstm1', 'sra1', 'srgap1', 'srgap3', 'srp72', 'srprb', 'srrm2', 'sry', 'ss18l1', 'sspop', 'ssr1', 'sstr3', 'sstr4', 'sstr5', 'ssx6p', 'ssx9p', 'st13p4', 'st14', 'st18', 'st20', 'st20-as1', 'st7l', 'stab2', 'stard3', 'stard5', 'stard7', 'stard9', 'stat3', 'steap1b', 'sth', 'sting1', 'stk16', 'stk17b', 'stk19', 'stk36', 'stra6', 'strada', 'strbp', 'strn4', 'stx11', 'stx18', 'stx1b', 'stxbp5', 'stxbp5l', 'styxl2', 'sucla2', 'sumo3', 'sumo4', 'sv2c', 'svep1', 'svil', 'sycn', 'syde2', 'syk', 'syn3', 'syngr4', 'synpo2', 'syt8', 'sytl3', 'tab2', 'tacc2', 'tacc3', 'tacr2', 'taf1l', 'taf5', 'taf7l', 'taf9', 'tanc1', 'tanc2', 'tapbp', 'tarbp2', 'tars1', 'tas1r2', 'tas1r3', 'tas2r1', 'tas2r10', 'tas2r13', 'tas2r4', 'tas2r40', 'tas2r41', 'tas2r5', 'tas2r7', 'tas2r8', 'tas2r9', 'tat', 'taz', 'tbx3', 'tbx6', 'tbxa2r', 'tbxas1', 'tcl1a', 'tcl1b', 'tdg', 'tdp2', 'tdrd1', 'tead3', 'tec', 'tek', 'tenm3', 'tenm4', 'tent5a', 'terf2', 'tert', 'tet3', 'tex11', 'tex14', 'tex15', 'tfb1m', 'tfrc', 'tgfa', 'tgfb1', 'tgfb1i1', 'tgfbi', 'tgfbr1', 'tgfbr2', 'thbd', 'thbs4', 'theg', 'thoc2', 'thoc3', 'thra', 'thsd1', 'tinf2', 'tjp3', 'tk2', 'tkt', 'tll1', 'tln2', 'tlr1', 'tlr10', 'tlr2', 'tlr3', 'tlr6', 'tlr9', 'tmc1', 'tmc2', 'tmc4', 'tmc6', 'tmc8', 'tmem132e-d', 'tmem14dp', 'tmem223', 'tmem230', 'tmem41b', 'tmem50a', 'tmem63c', 'tmem99', 'tmod2', 'tmod4', 'tmprss2', 'tmprss3', 'tnf', 'tnfaip3', 'tnfrsf10c', 'tnfrsf10d', 'tnfrsf11a', 'tnfrsf11b', 'tnfrsf13b', 'tnfrsf13c', 'tnfrsf14', 'tnfrsf17', 'tnfrsf4', 'tnfrsf8', 'tnfrsf9', 'tnfsf10', 'tnfsf11', 'tnfsf14', 'tnnc1', 'tnni3', 'tnrc18', 'tnrc6b', 'tob1', 'top1', 'top2a', 'top2b', 'top3a', 'top6bl', 'tox2', 'tp53bp1', 'tp53inp1', 'tpd52l1', 'tpmt', 'tpp1', 'tpsd1', 'tpsg1', 'tpt1', 'tpx2', 'tra', 'traf3ip2', 'traf3ip3', 'trappc14', 'trb', 'trbv7-9', 'trdn', 'treml4', 'trhr', 'trim10', 'trim2', 'trim31', 'trim38', 'trim40', 'trim47', 'trim49b', 'trim69', 'triobp', 'trip11', 'trmt6', 'trnt1', 'tsc1', 'tsen15', 'tspan17', 'tspan4', 'tspan6', 'tspan7', 'tspyl1', 'tspyl4', 'tssk2', 'tssk3', 'tst', 'ttbk1', 'ttc3', 'ttf2', 'ttk', 'ttll13p', 'ttmp', 'ttyh2', 'tubgcp2', 'tulp1', 'twf2', 'twist2', 'tyk2', 'tyro3', 'tyrobp', 'uaca', 'ubap1', 'ubd', 'ucp1', 'ufd1', 'ugt1a1', 'ulbp1', 'ulk1', 'umod', 'umodl1', 'umodl1-as1', 'unc119', 'unc5d', 'upb1', 'upf2', 'upf3a', 'upk1a', 'upk1b', 'urad', 'usp10', 'usp16', 'usp18', 'usp20', 'usp21', 'usp24', 'usp29', 'usp3', 'usp34', 'usp41', 'utf1', 'utrn', 'vamp1', 'vamp2', 'vangl1', 'vapa', 'vapb', 'vars1', 'vasp', 'vav1', 'vav2', 'vax2', 'vbp1', 'vcl', 'vcx3b', 'vdr', 'vit', 'vnn1', 'vnn2', 'vnn3', 'vpreb1', 'vps18', 'vps33a', 'vps35', 'vps41', 'vsx1', 'vwa3a', 'wars1', 'wasf3', 'wdr1', 'wdr11', 'wdr13', 'wdr17', 'wdr18', 'wdr19', 'wdr4', 'wfdc1', 'wfdc10b', 'wfdc3', 'wfdc5', 'wfdc8', 'wfdc9', 'wfs1', 'whamm', 'wipf3', 'wnk4', 'wnt1', 'wnt4', 'wrap73', 'wrn', 'wt1-as', 'wwox', 'xp32', 'yars1', 'ydjc', 'ywhah-as1', 'zbtb22', 'zbtb38', 'zbtb42', 'zc3h4', 'zcchc4', 'zfp36', 'zfpm1', 'zfpm2', 'zfr', 'zic3', 'zim2', 'zim3', 'zkscan8', 'zmat1', 'zmpste24', 'znf137p', 'znf154', 'znf200', 'znf207', 'znf211', 'znf214', 'znf215', 'znf219', 'znf22-as1', 'znf221', 'znf222', 'znf223', 'znf225', 'znf230', 'znf234', 'znf235', 'znf236', 'znf239', 'znf25', 'znf252p-as', 'znf256', 'znf263', 'znf264', 'znf266', 'znf273', 'znf287', 'znf304', 'znf334', 'znf337', 'znf354c', 'znf358', 'znf366', 'znf398', 'znf439', 'znf469', 'znf767p', 'znf806', 'znf818p', 'znf840p', 'znf844', 'zp2', 'zp3', 'zp4', 'zscan5c', 'zwint', '-']

for gene in missing_genes:
  download_cds(gene, '..\\New_CCDS')

#gene_dict = generate_gene_dict('..\\Table_S1_exonic.csv')
#chromos = get_chromo_seqs('..\\Full_Genome')
