#!/bin/bash

start_chr=22
stop_chr=22
# mix_pop="ASW"
pure_pop1="CEU"
pure_pop2="YRI"
logf="run_pipeline.log"
errf="run_pipeline.err"

touch $errf $logf
cat $errf >> $errf.full
cat $logf >> $logf.full
rm $errf $logf
# num_admixed_array=(10 200 200 250 250 300 300 350 350 400 400)
# num_pure_array=(300 30  20  30  20  30  20  30  20  30  20)

sample_ids_file="SampleIDs/igsr_samples.tsv"
num_admixed_array=(5)
num_pure_array=(200)
num_array_elems=${#num_admixed_array[@]}
num_recombinations=3
recombination_rate=0.001


## Main starts here
echo "$(date): Starting pipeline" | tee -a $logf >> $errf
bin/bashScripts/get_vcf_data.sh $start_chr $stop_chr>> $logf 2>> $errf
echo "$(date): get_vcf_data.sh done" | tee -a $logf >> $errf


# -----------------------

if [ ! -f $sample_ids_file ]; then
	echo "Missing ${sample_ids_file}"
	echo "Need to download the SampleIDs for all populations"
	exit 1
fi

# python bin/pythonScripts/write_sample_IDs.py ${sample_ids_file} ${mix_pop} ${pure_pop1} ${pure_pop2} >> $logf 2>> $errf
# echo "$(date): write_sample_IDs.py done" | tee -a $logf >> $errf

python bin/pythonScripts/write_sample_IDs.py ${sample_ids_file} ${pure_pop1} ${pure_pop2} >> $logf 2>> $errf
echo "$(date): write_sample_IDs.py done" | tee -a $logf >> $errf

# ----------------------

# bin/pythonScripts/select_populations.py ${start_chr} ${stop_chr} ${mix_pop} ${pure_pop1} ${pure_pop2}  >> $logf 2>> $errf
# echo "$(date): select_populations.py done" | tee -a $logf >> $errf

bin/pythonScripts/select_populations.py ${start_chr} ${stop_chr} ${pure_pop1} ${pure_pop2}  >> $logf 2>> $errf
echo "$(date): select_populations.py done" | tee -a $logf >> $errf


#@ToDo move below for STRUCTUREpainter part
#python bin/pythonScripts/write_sample_IDs.py ${sample_ids_file} ${pure_pop1} ${pure_pop2} >> $logf 2>> $errf
#bin/pythonScripts/select_populations.py ${start_chr} ${stop_chr} ${pure_pop1} ${pure_pop2}  >> $logf 2>> $errf

# bin/bashScripts/clean_chr_data_for_local_ancestry_split_new.sh ${start_chr} ${stop_chr} ${mix_pop} ${pure_pop1} ${pure_pop2} >> $logf 2>> $errf
bin/bashScripts/clean_chr_data_for_local_ancestry_split_new.sh ${start_chr} ${stop_chr} ${pure_pop1} ${pure_pop2} >> $logf 2>> $errf
echo "$(date): clean_chr_data_for_local_ancestry_split_new.sh done" | tee -a $logf >> $errf


mkdir -p pops_data/admixed > /dev/null 2> /dev/null
mkdir -p pops_data/admixture/bed > /dev/null 2> /dev/null

for (( index=0; index<$num_array_elems; index++ ))
do
	num_admixed=${num_admixed_array[$index]}
	num_pure=${num_pure_array[$index]}

	python bin/pythonScripts/create_admixed_chromosomes.py \
		--chr ${start_chr} \
		--num_admixed $num_admixed \
		--num_anchor $num_pure \
		--pops ${pure_pop1} ${pure_pop2} \
		--num_recombinations ${num_recombinations} >> $logf 2>> $errf
	
	# These lines will be run by run_plink_and_admixture.sh if running scripts individually
    admixed_output_homologous_file="pops_data/admixed/${pure_pop1}_${pure_pop2}_admixed_${num_admixed}admixed_${num_pure}pure_HOMOLOGOUS.vcf"
    admixed_output_allele_file="pops_data/admixed/${pure_pop1}_${pure_pop2}_admixed_${num_admixed}admixed_${num_pure}pure_ALLELE_vcf.txt"
    admixed_proportions_file="pops_data/admixed/${pure_pop1}_${pure_pop2}_admixed_${num_admixed}admixed_${num_pure}pure_proportions.txt"
    
    plink_output_file_prefix="pops_data/admixture/bed/${pure_pop1}_${pure_pop2}_admixed_${num_admixed}admixed_${num_pure}pure"
    # First we prune the file to eliminate SNPs that are linked
	plink2 --vcf ${admixed_output_homologous_file}  --indep-pairwise 50 5 0.5 --out tmp_plink1_out >> $logf 2>> $errf
	# Now lets use the prune.in file to generate the bed file
	plink2 --vcf ${admixed_output_homologous_file}  --extract tmp_plink1_out.prune.in --make-bed --thin 0.9 --out ${plink_output_file_prefix} >> $logf 2>> $errf
	/bin/rm tmp_plink1_out*

	admixture --cv ${plink_output_file_prefix}.bed 2 >> $logf 2>> $errf
	mv ${pure_pop1}_${pure_pop2}_admixed_${num_admixed}admixed_${num_pure}pure.2.* pops_data/admixture	

	#./create_all_admixed_chromosomes.sh >> $logf 2>> $errf
	echo "$(date): create_admixed_chromosomes.py done" | tee -a $logf >> $errf
	# ------------------------------

	python bin/pythonScripts/generate_test_set.py \
		--num_admixed ${num_admixed} \
		--chr ${start_chr} \
		--pops ${pure_pop1} ${pure_pop2} \
		--num_recombinations ${num_recombinations} >> $logf 2>> $errf

	# python bin/pythonScripts/generate_test_set.py >> $logf 2>> $errf
	echo "$(date): generate_test_set.py done" | tee -a $logf >> $errf
	# ------------------------------

	prefix="${pure_pop1}_${pure_pop2}_admixed_${num_admixed}admixed_${num_pure}pure"
	
	num_test_ids=-1	# -1 means all.
	python3 bin/pythonScripts/Test_Local_Ancestry_Inference_ASW_from_3Pops.py \
		--reference_filename pops_data/admixed/${prefix}_ALLELE_vcf.txt \
		--all_admix_filename pops_data/admixture/${prefix}.2.Q \
		--chrom_admix_filename pops_data/admixture/${prefix}.2.Q \
		--num_test $num_test_ids \
		--recombination_rate $recombination_rate \
		--test_filename test_input/${pure_pop1}_${pure_pop2}_${num_admixed}test_SNPs_ALLELE_vcf.txt >> $logf 2>> $errf
	
	#./test_ancestry.sh >> $logf 2>> $errf
	echo "$(date): test_ancestry.sh done" | tee -a $logf >> $errf
	# ------------------------------
	
	python3 bin/pythonScripts/evaluate_inferences.py \
		--inferences_filename results/source=${prefix}_ALLELE_vcf.txt_test=${pure_pop1}_${pure_pop2}_${num_admixed}test_SNPs_ALLELE_vcf.txt.txt \
		--true_ancestry_filename test_input/${pure_pop1}_${pure_pop2}_${num_admixed}true_population.csv \
		--save_plot_filename results/${prefix}_plot_5chroms_try2point5.png >> $logf 2>> $errf
		
	# ./evaluate_inferences.sh >> $logf 2>> $errf
	echo "$(date): evaluate_inferences.sh done" | tee -a $logf >> $errf
done
