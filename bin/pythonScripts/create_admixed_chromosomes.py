#!/usr/bin/env python3

'''
This script currently expects the output folder pops_data/admixed to be present.

This script should be run as part of create_all_admixed_chromsomes.sh. The arguments are there.

parses the output of `split_homologous_chr.py` and creates a number of admixed genomes from the two given populations
- run the script:
    mkdir -p pops_data/admixed
    python create_admixed_chromosomes.py --num_admixed_chromosomes 10 --num_anchor_chromosomes 300 --source_pops CEU YRI

Note: The above is almost what create_all_admixed_chromsomes.sh does, except that ALSO runs admixture and creates an
admixture folder in pops_data. So this should not really be run directly.

(outputs two files in data/admixed/ directory: a .vcf file containing admixed chromosomes, and a proportions.txt file
    containing the ancestry proportions of each of the generated chromosome, for comparison with the output
    of admixture program)

next steps:
- use plink to generate .bed file:
    mkdir -p data/admixture/new/bed
    plink --vcf data/admixed/CEU_YRI_admixed_n\=200.vcf --make-bed --out data/admixture/new/bed/CEU_YRI_admixed_n\=200

- run admixture on .bed file:
    admixture --cv data/admixture/new/bed/CEU_YRI_admixed_n\=100.bed 2
    mkdir -p data/admixture/new/admix
    mv CEU_YRI_admixed_n\=100.2.* data/admixture/new/admix
'''
import sys
import numpy as np
import pandas

import argparse
from collections import defaultdict
#@ToDo Commented the two imports below as they are unused
#    import itertools
#    import pdb

# Whats exactly is happening in the previous function? Lets rehash. We are passing two individuals from the two distinct populations to the
# function along with the Sample_Ids for each. In addition, we pass the entire allele.vcf data
# In the function, we first get randomly NUM_RECOMBINATION POS entries from the allele.vcf array in recombination_spots list
#
# Next, we identify the exact line number/position each recombination_spots POS entry is located in the allele.vcf file (see explanation
# before the previous LAMBDA statement).
# So, we have a list (chr_position_list) which has NUM_RECOMBINATION entries, with each entry corresponding to a line number in the allele.vcf file where
# the POS entry matches the entry in recombination_spots list.
# Next, we simply create two lists from this - a "start" list with has a 0 followed by the entries from "stop" and a new stop list
# with the entries from stop followed by the last line number of the allele file.
#
# We use the two list to iterate through the allele.vcf file to grab a section from start[i] to stop[i] , for the column corresponding to 
# the distinct individual tuple passed to the function, and saves it in the test_chr list. We alternate when picking the sections from the allele file
# from pop1 first, and then pop2, then pop1 again and so on. This results in a random "mixed indivdual SNPs " in the test_chr list. true_chr is either 0
# or 1 depending on whether the section is populated for pop1 or pop2. This true_chr is simply used to calculate the "proportion" of pop2 present in
# this "mixed individual" chromosome.
#
# So the function returns the test_chr which represents an admixed chromosome for the "mixed individual" created in the function from the two
# individuals from the two distinct populations passed to the function. In addition it also returns the ancestry proportion calculated.
#

#@ToDo Can we replace this with a logic such as
# 1. Get num_recombinations random numbers from 0 to len(chr_strands['POS']
# 2. For each random position from chr_strands['POS'], write the 'POS' value into recombination_spots tuple
# 3.       and write this random position as the 'index' (+1 for 1-index reference) in the stop variable (type = pandas.Series) below
def create_test_chromosome(chr_strands, pops, ref_IDs):
    ### sample from SNPs
    recombination_spots = (
        # This will grab num_recombinations unique values from the POS column in the allele.vcf file, which is read into
        # chr_strands using the pandas.parse_csv functions. replace=False ensures you do not get duplicates.
        # The list is sorted and re-indexed [0 .... num_recombinations] following the random sampling and sorting
        chr_strands['POS'].sample(num_recombinations, replace=False).sort_values().reset_index(drop=True)
    )

    # The command below is creating a pandas.Series object, with the entries as the 1-reference position numbers of the
    # ids from the recombination_spots list, in the chr_strands list.
    # For example, if recombination_spots[0] = 34403487 , then stop[0] would contain (index=0, value=53374), where
    # 53374 is the 1-indexed position (not 0 indexed) where 34403487 is found in the chr_strands array
    # The lambda function gets the list index in chr_strands where the condition is TRUE
    # The result of LAMBDA is (don't try to decode this ..... but this is what it does) that stop is a panda.series object
    #    with "num_recombinations" entries, where each entry is the position in the allele.vcf file or chr_strands list
    #    where the recombination_spots POS value is found in the file.
    # For e.g. if recombination_spots[0] == 34403487  then stop[0] = 53374
    # & if recombination_spots[1] == 36891858  then stop[1] = 62087
    # @ToDo Can we devise a simpler method to do this ..... this pandas --- lambda combo is very very confusing
    chr_position_list = pandas.Series(range(num_recombinations)).map(lambda x: np.where(chr_strands['POS'] == recombination_spots[x])[0][0] + 1)
    # The line below is simply adding a (0,0) entry at the start of the "stop" Series above
    start = pandas.concat([pandas.Series(0), chr_position_list], ignore_index=True)
    # The line below is simply adding a (len(chr_strands),len(chr_strands)) entry at the end of the "stop" Series above (not the start series)
    # Remember len(DataFrame) does not include the header line
    stop = pandas.concat([chr_position_list, pandas.Series(len(chr_strands))], ignore_index=True)
    # So at this point, both start and stop have the same number of entries

    # The next two lines are doing the example same thing. Creating a list of size len(chr_strands) with each value as 0
    # Two methods are being used to do the same thing - I guess somebody was learning and trying to do two things
    test_chr = np.zeros(len(chr_strands), dtype=np.uint8)
    true_chr = np.array([0] * len(chr_strands), dtype=np.uint8)
    for j in range(2):
        for i in range(j, len(start), 2):
            test_chr[start[i] : stop[i]] = chr_strands[ref_IDs[j]][start[i] : stop[i]]
            true_chr[start[i] : stop[i]] = j

    ancestry_proportions = true_chr.sum() / true_chr.shape[0]
    return test_chr, ancestry_proportions # proportion of SNPs selected from second ancestry


# Simple wrapper to call create_pure_chromosomes_from_one_pop function twice, once for each pure population list
# Returns a list with two arguments, first being a list of lists as pure chromosomes for a set of individuals from the two populations
# and 2nd as the list of individuals. We would always return two rows from this function - pop1 row and pop2 row.
def create_pure_chromosomes_from_all_pops(chr_strands, input_num_anchor_chromosomes, pop1_ids, pop2_ids):
    return [create_pure_chromosomes_from_one_pop(chr_strands, int(input_num_anchor_chromosomes/2.0), popIDs) for popIDs in [pop1_ids, pop2_ids]]


# This function is called for each set of distinct populations at a time. 
# The function will create a list of "pure" choromosome list for a random list of individuals from the population passed to it
# The assert statement below ensures that we use the minimum of the two lengths from the two populations, so that we get
# exactly num_anchor_chromosomes entries after two iterations of this function.
def create_pure_chromosomes_from_one_pop(chr_strands, half_input_num_anchor_chromosomes, popIDs):
    assert half_input_num_anchor_chromosomes <= len(popIDs)

    IDs_to_use = np.random.choice(popIDs, size=half_input_num_anchor_chromosomes, replace=False)
    chroms = [chr_strands[popID] for popID in IDs_to_use]
    return [chroms, IDs_to_use]


# This function is very similar to create_test_chromosome, except it returns three items
# 1st item is a list of lists - or list of "admixed chromosomes" - one for each distinct individual pair
#      where the "num_chromosomes" distinct individual pairs are randomly created in the function
# 2nd items is a list of ancestry proportions as a tuple - for pop1 and pop2
# and 3rd item is the list of random distinct individual pairs used in the function
# Please read create_test_chromosome to really understand what this function does, as this function
# is basically a "list" return of the result of that function.
def create_test_chromosome_set(chr_strands, pops_dict, num_chromosomes):
    ### pops_dict should map population name to list of haploid IDs corresponding to that population

    ID_pairs = []
    ### copy dict values (lists) since they will be modified; two copies of each haplotype, so ignore one
    ### remove IDs not in chr_strands
    pops_ordered = ['pop1', 'pop2']
    for i in range(num_chromosomes):
        if i % 20 == 0:
            print('on chromosome {}'.format(i))
        ID_pairs.append([pops_dict[x][np.random.randint(0, len(pops_dict[x]))] for x in pops_ordered]) ### get a random ID from each pop

    # The line below will populate the variable with a list of paired entries; first, the admixed chromosome list,
    # and second, the ancestry proportion in the list 
    test_chromosomes_and_proportions = [create_test_chromosome(chr_strands, pops_ordered, pair) for pair in ID_pairs]
    
    # The zip command below just extracts the two items from the list of pairs above into a list of admixed chromosomes
    # and the list of ancestry proportions
    test_chromosomes, pop2_ancestry_proportions = zip(*test_chromosomes_and_proportions)
    ancestry_proportions = [[1 - proportion, proportion] for proportion in pop2_ancestry_proportions]

    return test_chromosomes, ancestry_proportions, ID_pairs


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
    # Uncommenting the below will standardize the output with the same sequence of random numbers
    # np.random.seed(0)

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--chr', help='chromosome to start on.', type=int, default=22)
    parser.add_argument('--num_admixed', help='number of admixed chromosomes to create.', type=int)
    parser.add_argument('--num_anchor', help='number of pure chromosomes to include', type=int)
    parser.add_argument('--pops',
        help='population codes (e.g. "CEU YRI") for the two populations that will be used in creating the simulated test chromosomes',
        nargs=2, required=True)
    parser.add_argument('--num_recombinations', help='number of recombination events for each admixed chromosome', type=int)

    args = parser.parse_args()
    chr_number = args.chr
    num_admixed_chromosomes = args.num_admixed
    num_anchor_chromosomes = args.num_anchor
    pop1, pop2 = args.pops
    num_recombinations = args.num_recombinations

    populations = pop1 + "_" + pop2

    # homologous_filename = '/home/greg/School/popgen/data/chr22.phase3.{}.SNPs.homologous.txt'.format('ASW_CEU_YRI') ### TODO: change this

    #@ToDo the chr22.phase3 .... file should be renamed to just chr.phase3.... , because this is in the Chr21, Chr22 folder and so on
    allele_filename = 'pops_data/{}_Data/Chr{}/tmp/chr{}.phase3.{}.SNPs.allele.vcf'.format(populations, chr_number, chr_number, populations)

    pop1_ids_filename = 'SampleIDs/{}_Sample_IDs_haploid.txt'.format(pop1)
    pop2_ids_filename = 'SampleIDs/{}_Sample_IDs_haploid.txt'.format(pop2)

    # @ToDo might want to change this so the files don't go into the admixed folder in one big mess
    out_basename = 'pops_data/admixed/{}_admixed_{}admixed_{}pure'.format(populations, num_admixed_chromosomes, num_anchor_chromosomes)
    out_filename_homologous = out_basename + '_HOMOLOGOUS.vcf' ### 0|0 0|0 1|1 0|0
    out_filename_allele = out_basename + '_ALLELE_vcf.txt' ### 0 0 1 0

    proportions_out_filename = 'pops_data/admixed/{}_admixed_{}admixed_{}pure_proportions.txt'.format(populations, num_admixed_chromosomes, num_anchor_chromosomes)

    # Get the population id for each population (CEU and YRI) into two list variables
    pop1_ids = get_pop_IDs(pop1_ids_filename)
    pop2_ids = get_pop_IDs(pop2_ids_filename)

    #all_chr_strands = pandas.read_csv(homologous_filename, sep='\t', header=0, comment='#')
    #@ToDo the allele_filename is from the tmp/ folder, with no "hash" on the CHROM line. Can we instead use the file with the hash in the parent folder?
    all_chr_strands = pandas.read_csv(allele_filename, sep='\t', header=0, comment='#')
    all_chr_strands.drop_duplicates(inplace=True)

    pop1_ids, pop2_ids = map(lambda ids: list(filter(lambda ID: ID in all_chr_strands.columns, ids)), [pop1_ids, pop2_ids])

    # The logic below looks through the IDs in the population sample ids file, and finds the ones for which we have a
    # column in the allele file. This is because we may NOT have phase3 data for ALL the IDs for that population, but they
    # would all be listed in the Sample_IDs file for the population.
    # So, if we have 200 entries in CEU_Sample_IDs.txt but only 180 columns match, then our resultant ids would be 180.
    # 
    # The compicated lambda command is doing what the two (for ... for) loops are doing
    # resultant_pop1_ids = []
    # for id in pop1_ids:
    #     for chr_id in all_chr_strands.columns:
    #         if id == chr_id:
    #             resultant_pop1_ids.append(id)
    # pop1_ids = resultant_pop1_ids
    #
    #
    # resultant_pop2_ids = []
    # for id in pop2_ids:
    #     for chr_id in all_chr_strands.columns:
    #         if id == chr_id:
    #             resultant_pop2_ids.append(id)
    # pop2_ids = resultant_pop2_ids
    
    # So, at this point we have pop1_ids and pop2_ids that have corresponding columns in the allele & homologous .vcf files

    # The line below gets a list of lists, i.e. list of admixed chromosome list in the variable chromosomes
    # ancestry_proportions is a list of tuples with the proportion for each population
    # ID_pairs is a random list of indivual pairs for each of the two "distinct" proportions used to the get the other two variables in the function
    chromosomes, ancestry_proportions, ID_pairs = create_test_chromosome_set(all_chr_strands, {'pop1' : pop1_ids, 'pop2' : pop2_ids}, num_admixed_chromosomes)
    ### ancestry proportions is a list of [pop1_proportion, pop2_proportion] for each test chromosome
    print('created admixed chromosomes; creating anchor chromosomes')

    # The line below gets a list of pure chromosomes for the individuals from the two populations as two rows with
    # 1st arg of each row being a list of lists of pure chromosomes for a random list of pop ids, 
    # and 2nd arg being the list of population IDs randomly chosen in the function
    pure = create_pure_chromosomes_from_all_pops(all_chr_strands, num_anchor_chromosomes, pop1_ids, pop2_ids)
    # source1pure is chromosomes and pop ids from first population
    # source2pure is chromosomes and pop ids from second population
    source1pure, source2pure = pure

    # In the lines below, the ancestry prportions list obtained above for the admixed chromosomes, is extended to contain
    # num_anchor_chromosomes/2.0 entries with a tuple [1,0] implying this is a "population 1 pure proprotion" and the next line
    #extending to contain num_anchor_chromosomes/2.0 entries with a tuple [0,1] implying this a "population 2 pure proportion".
    num_anchor_chromosomes_per_ancestry = int(num_anchor_chromosomes/2.0)
    ancestry_proportions.extend([[1,0]] * num_anchor_chromosomes_per_ancestry)
    ancestry_proportions.extend([[0,1]] * num_anchor_chromosomes_per_ancestry)

    print('created anchor chromosomes; dropping unneeded columns')

    all_chr_strands.drop(all_chr_strands.columns[9:], axis=1, inplace=True)
    to_write = all_chr_strands
    to_write_allele = to_write.copy()

    print('finished dropping unneeded columns')

    print('building dataframe, starting with admixed chromosomes')

    # defaultdict is a neat way to create a dict where the "keys" will be dynamically created later
    # and when we check for a value for a key never previously allocated, it is automatically set to "0".
    used_names = defaultdict(int)
    # The statement below will assign 
    # index to 0, 1, 2,..... etc
    # chrom will have the admixed chromosome list 
    # and id_pair is the pair of individuals used to create the admixed chromosome list
    for index, (chrom, id_pair) in enumerate(zip(chromosomes, ID_pairs)):
        if index % 20 == 0:
            print('on column: {}'.format(index))

        # Using "name" to represent the name of the "offspring" from the two individuals
        name = '-'.join(id_pair).replace('_', '-')
        idx = used_names[name]      # idx will start with 0 due to the defaultdict use above
        used_names[name] += 1       # Record how many times we have a specific offspring name
        name = '{}-{}'.format(name, idx)    # Rename offspring name to include count to really make them unique

        # We removed all the columns after the 9th columns, so all the population id columns are dropped.
        # Now, we are adding the columns below for each admixed individual with the name as "HAxxxx-HAyyyy-n" (n >=0 )
        # and we are writing it as required for a .vcf file format (0|0  or 1|1 etc)
        to_write[name] = np.fromiter(('{0}|{0}'.format(snp) for snp in chrom), dtype='<U3')
        # Where homologous has the entries as "0|0" and allele has entries as "0" .
        # The line below simply adds the admixed chromosome list as a new column
        to_write_allele[name] = chrom

    if (num_anchor_chromosomes > 0):
        print('adding pure chromosomes to dataframe')
    else:
        print('not adding any pure chromosomes to dataframe')

    # The logic below is wanting to add all columns from source1pure chromosomes and source2pure chromosomes
    # The len(to_write.columns) and column variable is just so we can print which column we are currently writing
    for pure_chroms, IDs in [source1pure, source2pure]:
        for column, (chrom, ID) in enumerate(zip(pure_chroms, IDs), len(to_write.columns)):
            if column % 20 == 0:
                print('on column: {}'.format(column))

            to_write[ID] = np.fromiter(('{0}|{0}'.format(snp) for snp in chrom), dtype='<U3')
            to_write_allele[ID] = chrom

    print('finished building dataframe. Writing to .vcf file: {}'.format(out_filename_homologous))

    #with open(out_filename, 'w') as f:
    with open(out_filename_homologous, 'w') as f:
        f.write('##fileformat=VCFv4.2\n#')
        # f.write('##fileformat=VCFv4.2\n')
        # to_write.rename(columns={'CHROM' : '#CHROM'}, inplace=True)
        to_write.to_csv(f, sep='\t', index=False, mode='a')

    #with open(out_filename_homologous, 'w') as f_homologous:
    with open(out_filename_allele, 'w') as f_homologous:
        to_write_allele.to_csv(f_homologous, sep='\t', index=False)

    print('finished writing to .vcf file. writing ancestry proportions: {}'.format(proportions_out_filename))

    with open(proportions_out_filename, 'w') as f:
        if __name__ == '__main__':
            for ancestry_zero, ancestry_one in ancestry_proportions:
                f.write('{} {}\n'.format(ancestry_zero, ancestry_one))
    
    # At the end of this program, we should have a new homologous vcf file with columns for individuals with the admixed chromosomes, and pure 
    # chromosomes from population 1 and population 2. We will also have a similar "allele.txt" file with same number of columns
    # as the homologous.vcf file. 
    # We will also have a ancestry proportions file with the proportions for each individuals column
        