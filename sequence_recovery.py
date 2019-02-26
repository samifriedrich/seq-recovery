#!/usr/bin/env python

##################################################
# sequence_recovery
#
# Description:
# Compares alignments of chicken and zebra finch models to taeGut1, and
# creates a taeGut1 BED track of chicken sequence blocks that do not overlap
# anywhere in the corresponding zebra finch model
#
# Originally written in Python 3.6.7
##################################################
## Author: {Sami Friedrich}
## Copyright: Copyright {2019}, {seq_recovery}
## Credits: [{none}]
## License: {MIT}
## Version: {1}.{0}.{1}
## Maintainer: {Sami Friedrich}
## Email: {sami.r.friedrich@gmail.com}
##################################################


### Libs
import os
import subprocess
import pandas as pd
import re

# USER-DEFINED VARIABLES
# See documentation folder for these example files
my_dir = '~/Desktop/seq_recovery'  # directory where the following files can be found if not providing full paths
conversion_file = 'gene_to_transcript.txt'  # two-column text file linking ENS gene to transcript
ens_file = 'ENS_key.txt'  # text file with gene name, finch ENS gene(s), and chicken ENS gene(s)
psl_file = 'alignments-to-taeGut1.psl'  # psl file of chicken and finch alignments to taeGut1
bed_name = 'output_BED.txt'  # name of output BED track file, will output to my_dir unless otherwise specified

### Define standard psl headers
psl_headers = [
    'score',
    'mis-match',
    'rep-match',
    'Ns',
    'Q_gapcount',
    'Q_gapbases',
    'T-gapcount',
    'T-gapbases',
    'strand',
    'Qname',
    'Qsize',
    'Qstart',
    'Qend',
    'Tname',
    'Tsize',
    'Tstart',
    'Tend',
    'block_count',
    'blockSizes',
    'qStarts',
    'tStarts'
]

### Functions
# split but exclude empty elements
def split_no_empties(str, splitSym):
    return [s for s in str.split(splitSym) if s.strip() != '']

# range that includes the end position
def range1(start, end):
     return range(start, end+1)

# retrieve all psl data for a given gene
# NOTE: use requires that the 'ens_dict' dictionary is built and populated
def extract_psl_data(gene, model_key):

    models = ens_dict[gene][model_key]
    chromosome = []
    nBlocks = 0
    blockSizes = []
    tStarts = []
    tEnds = []

    # loop through each model's psl data and add to 'master' lists
    for model in models:
        psl_data = psl.loc[model]
        chromosome.append(psl_data['Tname'])
        nBlocks += psl_data['block_count']
        blockSizes.extend(split_no_empties(psl_data['blockSizes'], ','))
        tStarts.extend(split_no_empties(psl_data['tStarts'], ','))

    # check if chromosomes agree (relevant when there are multiple models)
    if chromosome[1:] != chromosome[:-1]:
        print(f'############ Chromosome incongruency for {gene} models {models}')
        print(chromosome)
    else:
        chromosome = chromosome[0]

    # convert coordinate and size values to ints
    blockSizes = list(map(int, blockSizes))
    tStarts = list(map(int, tStarts))

    # create list containing coordinates for where blocks end 'tEnds'
    for i in range(nBlocks):
        end = tStarts[i] + blockSizes[i]
        tEnds.append(end)

    return chromosome, nBlocks, blockSizes, tStarts, tEnds


### Build conversion table to translate between gene and transcript ID's

os.chdir(my_dir)
conversion_dict = {}

for line in open(conversion_file):
    line = line.strip('\n')
    line = line.replace(' ','')
    parts = line.split(',')
    gene_id = parts[0]
    conversion_dict.setdefault(gene_id, [])
    conversion_dict[gene_id] = parts[1]


### Build 'ens_dict' dictionary linking genes to their gene and transcript models in finch and chicken

ens_dict = {}


for line in open(ens_file):
    line = line.replace('"','')  # get rid of quotation marks
    line = line.replace('\n','')  # get rid of newline chars

    # only process genes with at least one ENS model in both chicken and zebra finch
    if line.count('ENSTG') < 1 or line.count('ENSGALG') < 1:
        #print('Missing a model in one or both species: ', line)
        continue

    parts = line.split('\t')

    # set key to gene name
    gene = parts[0]
    ens_dict.setdefault(gene, [])

    # get zf models and assign to key
    regex = r'\bENST\w+'
    zf_mod = re.findall(regex, line)

    # get chk models and assign to key
    regex = r'\bENSG\w+'
    chk_mod = re.findall(regex, line)

    ens_dict[gene] = {}
    ens_dict[gene]['zf_gene'] = zf_mod
    ens_dict[gene]['chk_gene'] = chk_mod

    # get zf transcript id's from conversion_dict
    trans_id = []

    for model in zf_mod:
        trans_id.append(conversion_dict[model])

    ens_dict[gene]['zf_trans'] = trans_id

    # get chicken transcript id's from conversion_dict
    trans_id = []

    for model in chk_mod:
        trans_id.append(conversion_dict[model])

    ens_dict[gene]['chk_trans'] = trans_id


### Remove header if present, read in psl data, and filter for top hits by query

for line in open(psl_file):
    # remove header
    if re.findall(r'^psLayout*', line):  # if first line begins with 'psLayout'
        # trim header (lines 1-5) off psl file:
        print(f'Trimming header from {psl_file}')
        cmd = "sed -e '1,5d' < " + psl_file + " > trimmed.psl"  # build command
        subprocess.check_output(cmd, shell=True)  # execute in shell
        psl_file = 'trimmed.psl'  # update psl_file to header-less trimmed.psl
        print(f'psl_file updated to: {psl_file}')

# read in data
psl_all = pd.read_csv(psl_file, sep='\t', header=None, names=psl_headers)
psl_all = psl_all.set_index('Qname')  # set index to query name

# filter for top-scoring alignments by query
max_filter = psl_all.groupby(['Qname'])['score'].transform(max) == psl_all['score']
psl = psl_all[max_filter]


### Block recovery

# initialize BED variables
chrom = []
chromStart = []
chromEnd = []
gene_name = []


for gene in ens_dict.keys():

    zf_chrom, zf_nBlocks, zf_blockSizes, zf_tStarts, zf_tEnds = extract_psl_data(gene, 'zf_trans')
    chk_chrom, chk_nBlocks, chk_blockSizes, chk_tStarts, chk_tEnds = extract_psl_data(gene, 'chk_trans')

    if zf_chrom != chk_chrom:
        print(f'############ Error: chromosomes do not agree between chicken and finch for {gene}')
        print(f'Zebra finch: {zf_chrom} \t Chicken: {chk_chrom}')

        # Exceptions:
        # TPCN3 split across chr3 (3') & chr3_random (5')
        # so overriding to chr3
        # RYR1 was excluded from analysis
        if gene != 'TPCN3':
            print('############ Excluded from exon recovery:', gene)
            continue # start over at next iteration of "for gene in ens_dict.keys()"

        if gene == 'TPCN3':
            zf_chrom = zf_chrom[1]  # set chromosome to chr3
            print('############ TPCN3 has been assigned zf_chrom[0]:', zf_chrom)

    nt_recov = 0

    for i in range(chk_nBlocks):

        add_to_BED = False

        chk_start = chk_tStarts[i]
        chk_end = chk_tStarts[i]+chk_blockSizes[i]
        chk_range = range1(chk_start, chk_end) # use custom function to include last bp in block (BED chromEnd is exclusive)
        chk_set = set(chk_range)

        for j in range(zf_nBlocks):

            zf_range = range1(zf_tStarts[j], zf_tStarts[j]+zf_blockSizes[j])

            if chk_set.intersection(zf_range):
                add_to_BED = False
                break
            else:
                add_to_BED = True

        if add_to_BED:
            chrom.append(zf_chrom)
            chromStart.append(chk_start)
            chromEnd.append(chk_end)
            gene_name.append(gene)
            nt_recov += chk_blockSizes[i]  # keep count of bases recovered

    # add total bases recovered to dictionary
    ens_dict[gene]['nt_recov'] = nt_recov


### Build dataframe and write to BED track

bed_track = pd.DataFrame(
    {'chrom': chrom,
     'chromStart': chromStart,
     'chromEnd': chromEnd,
     'name': gene_name
    })

bed_track.to_csv(bed_name, sep = ' ', header=False, index=False)
