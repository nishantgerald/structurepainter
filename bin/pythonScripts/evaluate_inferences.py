import argparse

import pandas as pd

import Chr_kmer_Informative_SNPs_Haploid_HMM as HMM

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--inferences_filename', help='path to VCF file containing inferences', required=True)
parser.add_argument('--true_ancestry_filename', help='path to VCF file containing true ancestry values', required=True)
parser.add_argument('--save_plot_filename', help='name to use for plot that will be generated (if not specified, plot will be shown instead of saved)', required=False, default=None)

if __name__ == '__main__':
	args = parser.parse_args()

	inferences_df = pd.read_csv(args.inferences_filename, sep='\t')
	true_df = pd.read_csv(args.true_ancestry_filename, sep='\t')

	ids = list(filter(lambda col: col != 'POS', inferences_df.columns))
	num_chromosomes = len(ids)

	print(HMM.evaluate(inferences_df, true_df))
	HMM.new_plot_ancestry_with_correct_results(inferences_df, true_df, image_filename=args.save_plot_filename)
