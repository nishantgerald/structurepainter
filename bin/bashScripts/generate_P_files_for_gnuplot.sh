# Generates the 4 P files to load into Gnuplot
# Should be run from Popgen folder
# @TODO make these generate for any populations, not just ASW/CEU/YRI

mkdir popgen_plots 2> /duv/null
# Finds the allele frequency difference between the populations in the reference chrosomomes after running ADMIXTURE
awk '{print $2 - $1}' pops_data/ASW_CEU_YRI_Data/Chr22/chr22.ASW_CEU_YRI.SNPs.homologous.25k.2.P | sed 's?^-??g' > popgen_plots/true_snps_diffs_chr22.P
# Remove any entries that are written in e notation (extremely small values because GNUPLOT can't deal)
grep -v '^[1-9]' popgen_plots/true_snps_diffs_chr22.P > popgen_plots/tmp_P
mv popgen_plots/tmp_P popgen_plots/true_snps_diffs_chr22.P
# Replace all differences less than 0.3 - making them appear as 0 in the plot, save to new file
sed "s?^0.[012].*?0?g" popgen_plots/true_snps_diffs_chr22.P > popgen_plots/true_snps_diffs_chr22_filtered.P

# Fine the allele frequency difference between the virtual admixed populations.
awk '{print $2 - $1}' pops_data/admixture/CEU_YRI_admixed_10admixed_300pure.2.P | sed 's?^-??g' > popgen_plots/CEU_YRI_admixed_diffs.P
# Remove any entries that are written in e notation (extremely small values because GNUPLOT can't deal)
grep -v '^[1-9]' popgen_plots/CEU_YRI_admixed_diffs.P > popgen_plots/tmp_P
mv popgen_plots/tmp_P popgen_plots/CEU_YRI_admixed_diffs.P
# Replace all differences less than 0.3 - making them appear as 0 in the plot, save to new file
sed "s?^0.[012].*?0?g" popgen_plots/CEU_YRI_admixed_diffs.P > popgen_plots/CEU_YRI_admixed_diffs_filtered.P