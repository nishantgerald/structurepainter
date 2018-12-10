### translates output from create_admixed_chromosomes into a format suitable for RFMix

### input: homologous VCF file, and ancestry proportions (each '1 0' or '0 1')
### output: 'alleles' and 'classes' files as described in RFMix input file documentation

import pandas as pd

import argparse
import os
import pdb

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('vcf_input_filename', type=str,
		help='name of file containing genetic data (with phased haplotypes separated into two separate columns, i.e. 0 1 instead of 0|1')

	parser.add_argument('ancestry_input_filename', type=str, \
		help='name of file containing ancestry proportions')

	args = parser.parse_args()

	alleles_filename = os.path.join('rfmix', 'alleles_{}'.format(os.path.basename(args.vcf_input_filename)))
	classes_filename = os.path.join('rfmix', 'classes_{}'.format(os.path.basename(args.vcf_input_filename)))

	vcf = pd.read_csv(args.vcf_input_filename, sep='\t')
	vcf.drop(vcf.columns[:9], axis=1, inplace=True) ### only retain haplotypes
	print('writing alleles to {}'.format(alleles_filename))
	with open(alleles_filename, 'w') as f:
		for _, series in vcf.iterrows():
			vals = series.map(str).str.cat()
			f.write('{}\n'.format(vals))

	print('writing classes to {}'.format(classes_filename))
	with open(args.ancestry_input_filename, 'r') as ancestry_input_f:
		with open(classes_filename, 'w') as classes_output_f:
			a1, a2, admixed = 0, 0, 0 ### keep track of how many of each type there were
			classes_list = []
			for line in ancestry_input_f:
				if line == '1 0\n':
					classes_list.append('1')
					a1 += 1
				elif line == '0 1\n':
					classes_list.append('2')
					a2 += 1
				else:
					classes_list.append('0')
					admixed += 1
			classes_output_f.write(' '.join(classes_list))

	print('found {} haplotypes from ancestry 1, {} from ancestry 2, and {} admixed haplotypes'.format(
		a1, a2, admixed))
