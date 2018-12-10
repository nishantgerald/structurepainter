#!/usr/bin/env python3

import pandas as pd

import pdb
from sys import argv
from os import path

if len(argv) <= 2:
    print("Usage: {0} <path to samples file> <population_1> [<population 2> .....]".format(argv[0]))
    exit()

input_samples_filename=argv[1]
pops = argv[2:]

ids_filepath = 'SampleIDs/{0}_IDs.txt'.format('_'.join(pops))

# Only make the IDs files if they don't already exist
if not path.isfile(ids_filepath):

    input_df = pd.read_csv(input_samples_filename, sep='\t')

    # Filter the dataframe to only match the pops given by user 
    # See learning_samples/pandas_dataframe.py to see how this works
    filtered_df = input_df[input_df['Population code'].isin(pops)]
    filtered_sample_name_df = filtered_df['Sample name']


    print('saving IDs belonging to any of ( {0} ) to {1}'.format(', '.join(pops), ids_filepath))
    # write IDs (from the Sample name column) from all the input populations in the ids_filepath defined above
    # index=False means the DataFrame index column is omitted. Default index is True
    filtered_sample_name_df.sort_values().to_csv(ids_filepath, index=False)

    # For each population, also make files for both chromsomes per individual (haploid files)
    for pop in pops:
        pop_df = input_df[input_df['Population code'] == pop]
        pop_sample_name_df = pop_df['Sample name']

        # For each Sample name, we append _1 and _2 to create the two haploids
        haploid_df = pd.concat([pop_sample_name_df.apply(lambda str: str + '_1'), pop_sample_name_df.apply(lambda str: str + '_2')]).sort_index()

        # Path the output filename for each population
        haploid_ids_filepath = 'SampleIDs/{0}_Sample_IDs_haploid.txt'.format(pop)
        print('saving haploid IDs for {0} to {1}'.format(pop,haploid_ids_filepath))

        # write haploid IDs (id_1, id_2) from one population
        haploid_df.sort_values().to_csv(haploid_ids_filepath, index=False)
else:
    print("{0} already exists".format(ids_filepath))