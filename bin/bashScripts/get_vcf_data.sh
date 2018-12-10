#!/usr/bin/env bash
# This script is expected to be run from the popgen folder.
# This script will download very large amounts of data and will take a LONG time to run.

URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
## Uncomment the line below if the above download link fails
# URL="http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3"
mkdir -p data/vcf
cd data/vcf

start_chr=$1
stop_chr=$2

for ((i=${start_chr};i<=${stop_chr};i++)); do
	vcf_filename="ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
	vcf_symlink="chr${i}.vcf.gz"

	# Download the vcf file for each chromosome in the data/vcf folder
	if [ ! -f ${vcf_filename} ]; then
		wget $URL/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	fi
	# Make a symlink with an easier name for each vcf file in the data/vcf folder
	if [ -f ${vcf_filename} -a ! -f ${vcf_symlink} ]; then
		ln -s $PWD/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz chr${i}.vcf.gz
	fi

done

