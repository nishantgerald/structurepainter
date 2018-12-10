from sys import argv
import time

# Called from: clean_chr_data_for_local_ancestry_split_new.sh
# This script must be executed from the popgen folder.
# This script takes roughly two minutes to run.

start_time = time.time()
print('starting split_homologous_chr')

in_filename = argv[1]
allele_filename = in_filename.replace('recode.vcf', 'allele.vcf')
homologous_filename = in_filename.replace('recode.vcf', 'homologous.vcf')

with open(in_filename, 'r') as in_file, open(allele_filename, 'w') as allele_file, open(homologous_filename, 'w') as homologous_file:
    header_found = False
    for line in in_file:
        if line[0] != '#': ### skip comment lines
            if not header_found: ### look for first line that isn't a comment
                header_found = True ### mark header as found
                line_list = line.split() ### header: NA001 NA002

                header_line = '\t'.join(line_list[0:9])
                # For each Sample_Id (NA001), make two entries (NA001_1 and NA001_2)
                for col in line_list[9:]:
                    header_line = header_line + '\t' + col + "_1"
                    header_line = header_line + '\t' + col + "_2"

                allele_file.write(header_line + '\n')
                homologous_file.write(header_line + '\n')
            else:
                line_list = line.split() ### Data line

                data_line = '\t'.join(line_list[0:9])
                # For each allele freq separated by |, we have two singlets, "0|1" is now "0 <tab> 1"
                for col in line_list[9:]:
                    allele_freq = col.split('|')    ### Split 0|1 into [0, 1]
                    data_line = data_line + '\t' + allele_freq[0]
                    data_line = data_line + '\t' + allele_freq[1]

                allele_file.write(data_line + '\n') ### 0 0 0 1

                data_line = '\t'.join(line_list[0:9])
                # For each allele freq separated by |, we have two pairs, "0|1" is now "0|0 <tab> 1|1"
                for col in line_list[9:]:
                    allele_freq = col.split('|')    ### Split 0|1 into [0, 1]
                    data_line = data_line + '\t' + allele_freq[0] + '|' + allele_freq[0]
                    data_line = data_line + '\t' + allele_freq[1] + '|' + allele_freq[1]

                homologous_file.write(data_line + '\n') ### 0|0 0|0 0|0 1|1

elapsed_time = time.time() - start_time
print('finished')
print('time elapsed: {}'.format(elapsed_time))
