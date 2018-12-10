#!/bin/bash

# This script must be executed from the popgen folder.
# Input: pops_data/ASW_CEU_YRI_Data/Chr22/chr22.phase3.ASW_CEU_YRI.SNPs.recode.vcf

# The script currently uses plink version 1.07. However, if you want to use plink 1.9, then uncomment the lines with
# plink2 below and ensure "plink2" points to the 1.9 version of plink.

# @ToDo: Why do we have the variable names as train, test1 and test2. Why not pop1, pop2 and pop3 ?
# @ToDo: Handle the linkage of pop input to the SampleIDs files so it is not always ASW_CEU_YRI
start_chr=$1
stop_chr=$2
shift 2
pops_list=$*


pops=$(echo $pops_list | sed "s/ /_/g")

pops_folder_base="pops_data/${pops}_Data"

for (( i=$start_chr; i<=$stop_chr; i++ )); do
	echo "###### Analyzing Chr$i ######"
	chr_folder_loc="${pops_folder_base}/Chr${i}"
	orig_chr_recoded_vcf="${chr_folder_loc}/chr${i}.phase3.${pops}.SNPs.recode.vcf"


	tmp_folder=${chr_folder_loc}/tmp
	if [ ! -d $tmp_folder ]; then
	    mkdir -p $tmp_folder
	fi

	working_chr_recoded_vcf="${tmp_folder}/chr${i}.phase3.${pops}.SNPs.recode.vcf"
	tmp_chr_header_file="${tmp_folder}/chr${i}.phase3.${pops}.header.txt"

    # The sed command below remove the hash from #CHROM line as the python script below expects it without the hash
	sed 's/#CHROM/CHROM/g' $orig_chr_recoded_vcf > $working_chr_recoded_vcf
	# Now, we extract ALL the (hash) header lines and save it in a working file, to be used later
	grep '#' $working_chr_recoded_vcf > $tmp_chr_header_file

	#### Split phased chromosomes
	# Output: pops_data/ASW_CEU_YRI_Data/Chr22/tmp/chr22.phase3.ASW_CEU_YRI.SNPs.recode.vcf
	python3 bin/pythonScripts/split_homologous_chr.py $working_chr_recoded_vcf

    # The two files below will be created in chr_folder_loc folder
	allele_filename=$(echo $working_chr_recoded_vcf | sed 's/recode.vcf/allele.vcf/g')
	homologous_filename=$(echo $working_chr_recoded_vcf | sed 's/recode.vcf/homologous.vcf/g')
	echo "finished running Python script"

	# What we have so far
	# - Input file - orig_chr_recoded_vcf - Chr22/chr22.phase3.......recode.vcf
	# - File without the "hash" - working_chr_recoded_vcf - Chr22/tmp/chr22.phase3......recode.vcf
	# - File that contains all the "hash" comment lines (except the "hash"CHROM line) - tmp_chr_header_file - Chr22/tmp/chr22.phase3....header.txt
	# - Allele file - allele_filename - Chr22/tmp/chr22.phase3.....allele.vcf
	# - Homologous file - homologous_filename - Chr22/tmp/chr22.phase3.....homologous.vcf

    # We now reinsert the header from the tmp header file we had saved earlier
	tmp_homologous_filepath_with_header="${tmp_folder}/chr${i}.${pops}.SNPs.homologous.withHeader.vcf"
	cat $tmp_chr_header_file $homologous_filename > ${tmp_homologous_filepath_with_header}

    # And, we add the (hash) in front of CHROM to rebuilt the vcf file after the changes made by the pythonscript above.
	homologous_filepath_with_header="${chr_folder_loc}/chr${i}.${pops}.SNPs.homologous.withHeader.vcf"
	sed 's/CHROM/#CHROM/g' ${tmp_homologous_filepath_with_header} > ${homologous_filepath_with_header}
	sed 's/CHROM/#CHROM/g' ${allele_filename} > ${chr_folder_loc}/chr${i}.phase3.${pops}.SNPs.allele.vcf
	echo "###### Finished splitting phased chromosomes ######"

    # @ToDo The vcftools command below is simply comverting the vcf file to a tped file format. Why do we need this if plink --vcf can accept a vcf file directly as input.
    # @ToDo Probably this is done because --vcf option is only available in plink 1.9 (not plink 1.07).
    # @ToDo Switch using plink 1.9 and eliminate the need for tped file generation
    vcftools_output_file="${chr_folder_loc}/chr${i}.${pops}.SNPs.homologous.withHeader2"
    vcftools --vcf ${homologous_filepath_with_header} --plink-tped --out $vcftools_output_file
	echo "finished vcftools"

	# This plink command creates a .prune.in file, which is used in the 2nd plink command below
	# The --indep-pairwise option has three values - window size in SNPs (e.g. 50), the number of SNPs to shift the window
    #       at each step (e.g. 5), the third parameter represents the r^2 threshold (where R^2 is the multiple correlation coefficient
    #       for a SNP being regressed on all other SNPs simultaneously - where 0 <= R^2 < 1)
    # So, a R^2 close to 1 results in higher threshold for linkage equilibrium, and thus lower number of SNPs that are pruned (i.e. larger prune.in file)
    plink_output_file_1="${tmp_folder}/chr${i}.${pops}.SNPs.homologous.LD"
    # plink2 --tfile ${vcftools_output_file} --indep-pairwise 50 5 0.5 --out ${plink_output_file_1}
    plink --noweb --tfile ${vcftools_output_file} --indep-pairwise 50 5 0.5 --out ${plink_output_file_1}
	echo "finished plink part1"

	# The plink_prune_file is created as an output of the above plink command.
	plink_prune_file="${tmp_folder}/chr${i}.${pops}.SNPs.homologous.LD.prune.in"
    plink_output_file_2="${chr_folder_loc}/chr${i}.${pops}.SNPs.homologous.25k"
    # plink2 --tfile ${vcftools_output_file} --extract ${plink_prune_file} --thin-count 25000 --make-bed --out ${plink_output_file_2}
    plink --noweb --tfile ${vcftools_output_file} --extract ${plink_prune_file} --thin .9 --make-bed --out ${plink_output_file_2}   # --thin 0.9 matches the total count to be close to --thin-count 25000
	echo "###### Selected 25k Independent SNPs ######"

    # plink_file_2_bed is generated by the previous plink command
    # Run the following on the command line to find the argument "K" value (shown as 2 below)
    # for k in 1 2 3 4 5
    # do
    #   admixture --cv chr22.ASW_CEU_YRI.SNPs.homologous.25k.bed $k | grep 'CV error'
    # done
    # Look for Output like .... CV error (K=2): 0.83889 .... and use the value of K where CV error is lowest.
    
    plink_file_2_bed="${plink_output_file_2}.bed"
    logf="${chr_folder_loc}/log.out"
    admixture --cv ${plink_file_2_bed} 2 | tee ${logf}
    mv chr${i}.${pops}.SNPs.homologous.25k.2.* ${chr_folder_loc}/
	echo "###### Completed ADMIXTURE estimates ######"

	# rm -rf ${tmp_folder}
done