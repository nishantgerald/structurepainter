'''
This file is similar to create_admixed_chromsomes.py, which has additional comments

Passes the output of s='split_homologous_chr.py', and creates a number of admixed genomes from the two given populations

This script generates two output tsv files in popgen/test_input/
    First output (true_population.csv) has true ancestry information (pop1, pop2). Refers to which true
     population the SNP was picked from to create the virtual admixed individual
    Second output (test_SNPs_ALLELE_vcf.txt) is a subset of the vcf file with the specified number of admixed individuals.
     Thus, it has 0 or 1 for each SNP (reference or alt) for each admixed individual. Each admixed invidividual was
     created from 2 homologous "parents" (really one of the two split chromsomes of each original parent individual)
     Whether the SNP was from a pop1 or pop2 individual is found in the first output file
'''

import numpy as np
import pandas as pd

import itertools
import pdb
import argparse
import os

#### Function to create test chromosomes ####
#### Create "admixed" populations from two pure populations
def create_test_chromosomes(chr_strands, pops, ref_IDs, num_recombinations):
    test_pos_chromosome = chr_strands['POS'].sample(num_recombinations, replace=False) ### sample from SNPs
    assert(len(test_pos_chromosome) == num_recombinations)

    test_pos_chromosome.sort_values(inplace=True)
    test_pos_chromosome.reset_index(drop=True, inplace=True)

    test_chr_idx_stop = pd.Series(range(num_recombinations)).map(lambda x: np.where(chr_strands['POS'] == test_pos_chromosome[x])[0][0] + 1)
    test_chr_idx_start = pd.concat([pd.Series(0), test_chr_idx_stop], ignore_index=True)
    test_chr_idx_stop = pd.concat([test_chr_idx_stop, pd.Series(len(chr_strands))], ignore_index=True)
    start = test_chr_idx_start
    stop = test_chr_idx_stop

    test_chr = np.zeros(len(chr_strands), dtype=int)
    true_chr = np.array([''] * len(chr_strands), dtype='<U4') ### need to specify max string length
    for j in range(len(pops)):    # pops = 2
        for i in range(j, len(test_chr_idx_start), len(pops)):    #
            test_chr[start[i] : stop[i]] = chr_strands[ref_IDs[j]][start[i] : stop[i]]
            true_chr[start[i] : stop[i]] = pops[j]

    return [test_chr, true_chr]


def create_test_chromosome_set(chr_strands, pops_dict, num_chromosomes, num_recombinations):
    ### pops_dictionary should map population name to list of haploid IDs corresponding to that population
    ### if num_chromosomes == -1, as many chromosomes as possible will be created

    ID_pairs = []
    ### copy dict values (lists) since they will be modified; two copies of each haplotype, so ignore one
    ### remove IDs not in chr_strands
    # pops_dict = {key: list(filter(lambda v: v[-2:] != '_2' and v in chr_strands.columns, val)) for key, val in pops_dictionary.items()}

    pops_ordered = ['pop1', 'pop2']
    for i in itertools.count():
        ### keep getting a random ID from each population, until either one of the ID lists is empty, or num_chromosomes is exceeded
        if [] in pops_dict.values() or (num_chromosomes != -1 and i >= num_chromosomes):
            break
        ID_pairs.append([pops_dict[pop].pop(np.random.randint(0, len(pops_dict[pop]))) for pop in pops_ordered]) ### get a random ID from each pop

    return [create_test_chromosomes(chr_strands, pops_ordered, pair, num_recombinations) for pair in ID_pairs], ID_pairs


def get_pop_IDs(filename):
    ret = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip('\n')
            if line != '':
                ret.append(line)
    return ret


if __name__ == '__main__':
    np.random.seed()
    # Uncommenting the below will standardize the output
    # np.random.seed(0)

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--num_admixed', help='number of admixed chromosomes.', type=int, default=10)
    parser.add_argument('--chr', help='chromosome to start on.', type=int, default=22)
    parser.add_argument('--pops',
        help='population codes (e.g. "CEU YRI") for the two populations that will be used in creating the simulated test chromosomes',
        nargs=2, required=True)
    parser.add_argument('--num_recombinations', help='number of recombination events for each admixed chromosome', type=int)


    args = parser.parse_args()
    num_admixed_chromosomes = args.num_admixed
    chr_number = args.chr
    sourcepop1, sourcepop2 = args.pops
    num_recombinations = args.num_recombinations

    populations = sourcepop1 + "_" + sourcepop2

    sourcepop1_ids_filename = 'SampleIDs/{}_Sample_IDs_haploid.txt'.format(sourcepop1)
    sourcepop2_ids_filename = 'SampleIDs/{}_Sample_IDs_haploid.txt'.format(sourcepop2) 
    
    #@ToDo We need to incorporate chr_number somewhere in our output folder path.
    # File below is from the tmp folder because the Chrom line (header line) does not and SHOULD not have a "HASH"
    allele_filename = 'pops_data/{}_Data/Chr{}/tmp/chr{}.phase3.{}.SNPs.allele.vcf'.format(populations, chr_number, chr_number, populations)

    sourcepop1_ids = get_pop_IDs(sourcepop1_ids_filename)
    sourcepop2_ids = get_pop_IDs(sourcepop2_ids_filename)

    all_chr_strands = pd.read_csv(allele_filename, sep='\t', header=0, comment='#')
    # source_df = pd.read_csv(allele_filename, sep='\t', header=0, comment='#')
    all_chr_strands.drop_duplicates(inplace=True)

    # The function below scans all_chr_strands (the input allele.vcf file across the populations of interest), and extracts the respective
    # pop ids from pop1 and pop2 in the two lists (as sourcepop1_ids and sourcepop2_ids)
    sourcepop1_ids, sourcepop2_ids = map(lambda ids: list(filter(lambda ID: ID in all_chr_strands.columns, ids)),[sourcepop1_ids, sourcepop2_ids])
    
    test_chr = create_test_chromosome_set(all_chr_strands, {'pop1' : sourcepop1_ids, 'pop2' : sourcepop2_ids}, num_admixed_chromosomes, num_recombinations)

    ancestry_df = pd.DataFrame(all_chr_strands[all_chr_strands.columns[:9]])
    # homologous_df = pd.DataFrame(all_chr_strands[all_chr_strands.columns[:9]])
    chrom_df = pd.DataFrame(all_chr_strands[all_chr_strands.columns[:9]])

    for (test, true), id_pair in zip(*test_chr):
        name = '-'.join(id_pair).replace('_', '-')
        ancestry_df[name] = true
        chrom_df[name] = test
        # homologous_df[name] = "{0}|{0}".format(test)

    # @Todo should we create a homolgoous vcf file as well for consistency? We only make an allele file here.

    # chrom_df.to_csv('test_input/CEU_YRI_test_SNPs_ALLELE_vcf.txt', sep='\t', index=False)
    if not os.path.exists('test_input'):
        os.makedirs('test_input')
    ancestry_df.to_csv('test_input/{}_{}_{}true_population.csv'.format(sourcepop1, sourcepop2, num_admixed_chromosomes), sep='\t', index=False)
    chrom_df.to_csv('test_input/{}_{}_{}test_SNPs_ALLELE_vcf.txt'.format(sourcepop1, sourcepop2, num_admixed_chromosomes), sep='\t', index=False)
    # homologous_df.to_csv('test_input/CEU_YRI_test_SNPs_homologous.vcf', sep='\t', index=False)
