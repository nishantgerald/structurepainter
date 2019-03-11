# STRUCTUREpainter
#### Software Needed ####

Need to install:
- [admixture](http://www.genetics.ucla.edu/software/admixture/download.html),
- [plink](http://zzz.bwh.harvard.edu/plink/download.shtml).
- [plink2](https://www.cog-genomics.org/plink2)
- [vcftools](https://vcftools.github.io/downloads.html).

Add these programs to your PATH.

#### Steps to Install on MacOS ####

- vcftools install
    - brew install vcftools
    
(if previous fails, compile from source using the steps below)
    - Install XCode Command Line Tools from developer.apple.com
    - brew install pkg-config
    - git clone https://github.com/vcftools/vcftools.git
    - cd vcftools
    - ./autogen.sh
    - ./configure
    - make
    - make install

- admixture install
    - Download from the admixture website, the MacOS binary and uncompress it
    - Make a symlink to admixture in your /usr/local/bin so it is in your $PATH
    - Admixture manual - https://www.genetics.ucla.edu/software/admixture/admixture-manual.pdf

- plink install
    - Download from the plink website, the MacOS binary and unzip it
    - Make a symlink to plink in your /usr/local/bin so it is in your $PATH

- plink2 install
    - Download from the plink2 website, the MacOS binary and unzip it
    - Make a symlink to plink2 in your /usr/local/bin so it is in your $PATH

#### Required Python packages: ####
- ggplot
- numpy
- pandas
- pymongo
- collections
- itertools
- pdb
- argparse
- matplotlib
- scipy
- time
- seaborn

#### Directory Structure ####

The following directories will be created when you download from git
- **bin/bashScripts** : This folder contains all the bash scripts used in this pipeline. The only bash script not in this folder is run_pipeline.sh, which is present in the top level folder.

- **bin/pythonScripts** : This folder contains all the python scripts used in the pipleline.

- **bin/RScripts** : This folder contains the RScripts, which have been converted to python. So, the RScripts are only here for reference.

- **docs** : This folder has some additional help documents related to the pipeline.
  
- **SampleIDs**: This folder will contain the file with the sample ids for the populations of interest. The file is already present when downloading from Git. However, you can download the file yourself from the IGSR website.
  
The following directories will be created during the pipeline processing
- **data/vcf** : This folder will contain the downloaded phased vcf file from
    - ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/. (Mirror: http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/)

    The vcf file can be downloaded using the script get_vcf_data.sh, or manually downloaded and move into data/vcf/ folder.

- **pops_data** : This folder will contain subfolders for populations of interest and the results from admixture processing in the pipeline.
    - As an example, you may see the following folders
        - pops_data/ASW_CEU_YRI_Data
        - pops_data/CEU_YRI_Data
        - pops_data/admixture
        - pops_data/admixed
        
- **results** : This folder will contain the output of StructurePainter when executed in the pipeline.

- **test_input** : This folder will contain the output from the generate_test_set.py program run in the pipeline.

### Running STRUCTUREPAINTER ###

There are currently two ways to run STRUCTUREPAINTER. 
1. You can run it as a single pipeline using run_pipeline.sh. 
2. Or you can run the scripts individually.

#### Preparing Input Data artifacts ####
The pipeline requires two input data artifacts. 
1. The vcf file with the phased chromosome details from the 1000 Genome project
2. And the sample ids for the populations from the IGSR (The International Genome Sample Resource) site (http://www.internationalgenome.org/data-portal/sample).

While the vcf file can be downloaded using a script provided here (bashScripts/get_vcf_data.sh), the sample IDs need to
be downloaded manually. The file in tsv (Tab Separated Value) format is obtained from IGSR (http://www.internationalgenome.org/data-portal/sample )
Simply download the unfiltered list from the website and save it in the folder SampleIDs. A previously downloaded version is present in Git.

#### Running STRUCTUREPAINTER as a pipeline ####

ALL scripts below must be run from the top level `popgen` folder

Running it with run_pipeline.sh is detailed below:

_@ToDo: Break up run_pipeline.sh into two, the 2nd named run_verification.sh, and clean the two scripts._

~~~
The pipeline script can be run using the command ./run_pipeline.sh.

However, you would need to edit the script and change the variables at the start
of the script for the particular requirements of the pipeline. These variables are

- start_chr AND stop_chr: These two variables will limit the pipeline to a subset of chromosome numbers
- pure_pop1 AND pure_pop2: These two variables define the 3 character designation for the two 'pure populations' of interest
- logf AND errf: These two variables have two files which will contain the log and error output from the pipeline process.
- sample_ids_file: This variable is the file location for the samples ids file downloaded from 1000 genomes.
                   (Previously downloaded version in Git - SampleIDs/igsr_samples.tsv)
- num_admixed_array: This array variable sets the number of admixed individuals to be generated for use by STRUCTUREPAINTER. 
                     It is to be used in conjunction with num_pure_array variable, and should have the same number of entries.
- num_pure_array: This array variable sets the number of pure individuals to be generated for use by STRUCTUREPAINTER. 
                  It is to be used in conjunction with num_admixed_array variable, and should have the same number of entries.
- num_recombinations: This variable is the number of recombination ("crossing-over") events used to create the virtual admixed chromosomes
- recombination_rate: This variable is used to modify the recombination rate used by STRUCTUREPAINTER.
                      @ToDo This should be taking a recombination map instead of one single rate

Once the variables above have been adjusted for the pipeline, ./run_pipeline.sh will perform all the
operations from downloading the vcf data (if not already downloaded) to running STRUCTUREPAINTER.

The final output from the pipeline will be present in the results directory.
~~~


#### Running STRUCTUREPAINTER with individual scripts ####

ALL scripts below must be run from the top level `popgen` folder

~~~
1. bin/bashScripts/get_vcf_data.sh
~~~

Download the vcf files for the chromosomes of interest
    
   This file will download the vcf files for the chromosomes (edit the file if not ALL chromosomes should be downloaded)
   and copies the files in the folder data/vcf. The script will also create a "shortened symbolic link" (e.g. chr22.vcf.gz), which is used in
   the select_populations.py script below.

~~~
2. python bin/pythonScripts/write_sample_IDs.py {sample_ids_file} {pure_pop1} {pure_pop2}
~~~

This script will use the unfiltered Sample IDs (from the SampleIDs folder) and extract the population IDs only for the population types specified as input arguments to the script.

   - `sample_ids_file`: This variable is the file location for the samples ids file downloaded from 1000 genomes.
                       (Previously downloaded version in Git - SampleIDs/igsr_samples.tsv)
   - `pure_pop1 AND pure_pop2`: These two variables define the 3 character designation for the two 'pure populations' of interest


~~~
3. bin/pythonScripts/select_populations.py ${start_chr} ${stop_chr} ${pure_pop1} ${pure_pop2}
~~~

Recodes the VCF file for each chromosome to a new VCF file of the pure populations.
    
   The script expects the following input files to be present.
   - `data/vcf/chr22.vcf.gz` - This is the vcf data file (one for each chromosome number). This file is created with the
                                proper name if the script `get_vcf_data.sh` is used to download the vcf file. Else,
                                the file should be renamed as shown here.
   - `SampleIDs/CEU_YRI_IDs.txt` - This is the file created by the write_sample_IDs.py script above, for the populations of interest.

   Script arguments:
   - `start_chr and stop_chr`: Start and stop chromosome number to be used for processing in the pipeline.
   - `pure_pop1 AND pure_pop2`: These two variables define the 3 character designation for the two 'pure populations' of interest

   Script runtime for Chr22: about 8 minutes.

   Calls: `vcftools`


~~~
4. bin/bashScripts/clean_chr_data_for_local_ancestry_split_new.sh ${start_chr} ${stop_chr} ${pure_pop1} ${pure_pop2}
~~~

Generate the homologous and allele files for the populations of interest from the outputs of the previous script. Also, run admixture to create the .Q files from the allele file (to be used later by STRUCTUREPAINTER, see Admixture documentation for details).
   

   The script expects the following input files to be present.
   - `pops_data/*_Data/Chr22/*.recode.vcf` - The recoded VCF file from select_populations.py is located in this folder for a specific chromosome.
   
   Script arguments:
   - `start_chr and stop_chr`: Start and stop chromosome number to be used for processing in the pipeline.
   - `pure_pop1 AND pure_pop2`: These two variables define the 3 character designation for the two 'pure populations' of interest

   Script runtime for Chr22: about 10 minutes.

   Calls: 
        `bin/pythonScripts/split_homologous_chr.py`
        `vcftools`, `plink`, `admixture` 
   
   The output from the script can be found in the folder pops_data/*_Data/Chr22/

   _@ToDo: The pops_data/*_Data/Chr*/tmp/ needs to be cleaned up. What files can be deleted or renamed? The tmp folder files are being used below at this time._


~~~
5. bin/pythonScripts/create_admixed_chromosomes.py --chr ${start_chr}  --pops ${pure_pop1} ${pure_pop2} --num_admixed ${num_admixed} --num_anchor ${num_pure}--num_recombinations ${num_recombinations}
~~~
Parses the output of the previous script (`clean_chr_data_for_local_ancestry_split_new.sh`, which calls `split_homologous_chr.py`) and creates an admixed training set. This is where we use the list of sampleIds per population created by write_sample_IDs.py.

   The script expects the following input files to be present.
   - Following folders must be present: pops_data/admixed  and   pops_data/admixture/bed (@ToDo: to create empty folders in git, make a .keep file in the folder.)
   - In `SampleIDs/` directory: two lists of sample IDs, one for each training population, created by the write_sample_IDs.py script above. e.g. SampleIDs/CEU_YRI_IDs.txt
   - In `pops_data/*_Data/Chr*/tmp/` directory: source of genetic data, e.g. `chr22.phase3.ASW_CEU_YRI.SNPs.allele.vcf`

   Script arguments:
   - `start_chr`: Start chromosome number to be used for processing in the pipeline.
   - `pure_pop1 AND pure_pop2`: These two variables define the 3 character designation for the two 'pure populations' of interest
   - `num_admixed`: This variable sets the number of admixed individuals to be generated for use by STRUCTUREPAINTER.
   - `num_pure`: This variable sets the number of pure individuals to be generated for use by STRUCTUREPAINTER
   - `num_recombinations`: This variable is the number of recombination ("crossing-over") events used to create the virtual admixed chromosomes

   Script runtime for Chr22: about 1 minute

   Output data:
   - In `pops_data/admixed/` directory:
     - admixed training set, in two forms:
       - e.g. `CEU_YRI_admixed_40admixed_100pure_ALLELE_vcf.txt`
       - e.g. `CEU_YRI_admixed_40admixed_100pure_HOMOLOGOUS.vcf`
         - These have the same column names, but when the allele file has `0`, the homologous vcf has `0|0`, and when the allele file has `1`, the homologous vcf has `1|1`.
       - True ancestry proportions for the admixed training set
         - e.g. `CEU_YRI_admixed_40admixed_100pure_proportions.txt`
         - **not** used by STRUCTUREpainter (we only have this information here because of our artificially created training set). 
         - This is included as verification. The advantage of STRUCTUREpainter is that pure populations are NOT needed.

~~~
6. bin/bashScripts/run_plink_and_admixture.sh ${pure_pop1} ${pure_pop2} ${num_admixed} ${num_pure}
~~~

Runs admixture on the admixed training set.

OPTIONAL: Use `plot_admix_results.py` script to see how well ADMIXTURE does in inferring ancestry proportions.
 
   Script arguments:
   - `pure_pop1 AND pure_pop2`: These two variables define the 3 character designation for the two 'pure populations' of interest
   - `num_admixed`: This variable sets the number of admixed individuals to be generated for use by STRUCTUREPAINTER.
   - `num_pure`: This variable sets the number of pure individuals to be generated for use by STRUCTUREPAINTER

   Script runtime for Chr22: about 1 minute

   Output data:
   - in `pops_data/admixture/` directory:
     - e.g. CEU_YRI_admixed_5admixed_200pure.2.P and CEU_YRI_admixed_5admixed_200pure.2.Q


~~~
7. bin/pythonScripts/generate_test_set.py --chr ${start_chr}  --pops ${pure_pop1} ${pure_pop2} --num_admixed ${num_admixed} --num_recombinations ${num_recombinations}
~~~

This script is similar to create_admixed_chromosomes.py. However, this script generates the virtual admixed populations as two tsv (tab separated values) files (see Output data below).

   The script expects the following input files to be present.
   - In `SampleIDs/` directory: two lists of sample IDs, one for each training population, created by the write_sample_IDs.py script above. e.g. SampleIDs/CEU_YRI_IDs.txt
   - In `pops_data/*_Data/Chr*/tmp/` directory: source of genetic data, e.g. `chr22.phase3.ASW_CEU_YRI.SNPs.allele.vcf`


   Script arguments:
   - `start_chr`: Start chromosome number to be used for processing in the pipeline.
   - `pure_pop1 AND pure_pop2`: These two variables define the 3 character designation for the two 'pure populations' of interest
   - `num_admixed`: This variable sets the number of admixed individuals to be generated for use by STRUCTUREPAINTER.
   - `num_recombinations`: This variable is the number of recombination ("crossing-over") events used to create the virtual admixed chromosomes

   Script runtime for Chr22: about 10 seconds

   - Output data:
     - In `test_input/` directory (test input to STRUCTUREpainter):
       - csv containing ancestry information e.g. CEU_YRI_5true_population.csv
       - csv containing genetic information e.g. CEU_YRI_5test_SNPs_ALLELE_vcf.txt

~~~
8. bin/pythonScripts/Test_Local_Ancestry_Inference_ASW_from_3Pops.py \
        --reference_filename ${admixed_allele_vcf} --all_admix_filename ${all_admixture_Q} --chrom_admix_filename ${chrom_admixture_Q} --test_filename {genetic_info_tsv_file} \
        --num_test ${num_test_ids} \
        --kmer ${window_size} --num_windows ${num_sliding_windows} --seed {random_seed}`
~~~
This is the STRUCTUREPainter script, which processes the data generated previously in the pipeline. Run the local ancestry method which selects ancestry informative SNPs, estimates the transition + emission matrices, and iterates through the Hidden Markov Model (HMM).

   Script arguments:
   
   Required Arguments:
   - `reference_filename` (admixed_allele_vcf) - path to file that contains training set - ALLELE_vcf.txt from pops_data/admixed folder.
   - `all_admix_filename` (all_admixture_Q) - path to overall admix filename - admixture output .Q file from pops_data/admixture folder. This needs to be a combined Q file for all 22 autosomal chromosomes.
   - `chrom_admix_filename` (chrom_admixture_Q) - path to admix filename for this chromosome - admixture output .Q from pops_data/admixture folder
   - `test_filename` (genetic_info_tsv_file) - path to file that contains chromosomes for which ancestry should be inferred - ALLELE.vxf.tx tsv from test_input folder
   
   Optional Arguments:
   - `num_test` (num_test_ids)- number of test chromosomes to evaluate (**default == -1 (which means all)**)
   - `kmer` (window_size) - Window size to use. This is also used as the number of sliding windows (**default 5**)
   - `num_windows` (num_sliding_windows) - number of sliding windows to use. If not specified, the value for kmer is used (**default 5**)           
      - Number of ways to split each chromosome into windows. If specified, should be less than or equal to `kmer`.
      - Another way to think about what this parameter means: for each SNP, how many calls do you want to make and then average together to determine ancestry? `num_windows` exactly determines this.
      - Using values other than `kmer` is mostly untested right now.
      - Note that approximately, runtime should scale linearly with each of `kmer` and `num_windows`.
   - `seed` (random_seed) - Any number to be used to seed the random number generator (**default 0**)
      - Set to `None` for irreproducible results (varying randomly from run to run), or to a different integer value to see a different set of results that will be the same from run to run.

   Script runtime for Chr22: about xx minutes

   - Output data:
       - The output file is created in the results folder, which will be used by evaluate_inferences.py to generate the ancestry plot. 
    
   _@ToDo Generate the input file arguments within the script using known information (chromsome number and populations) and known directory structure like previous scripts._
   
   _@ToDo The logic for passing a combined 22 chromsome VCF file to ADMIXTURE needs to happen__
    
~~~
9. evaluate_inferences.py --inferences_filename ${STRUCTUREPainter_output) --true_ancestry_filename (ancestry_information_csv) --save_plot_filename (output_png_file)
~~~
Will plot the results of STRUCTUREpainter against the true ancestry values for the test set of chromosomes.

   Script arguments:
   
   Required Arguments:    
   - `inferences_filename` (STRUCTUREPainter_output) - Output of StructurePainter from the results/ folder.
   - `true_ancestry_filename` (ancestry_information_csv) - csv containing true ancestry information from the test_input folder e.g. CEU_YRI_5true_population.csv.
   
   Optional Arguments:
   - `save_plot_filename` (output_png_file) - name to use for plot that will be generated (if not specified, plot will be shown instead of saved).


        ----Stopped here when editing Readme--------

#### Example running individual scripts ####
    - example usage `python ./bin/pythonScripts/write_sample_IDs.py SampleIDs/igsr_samples.tsv ASW CEU YRI`
        - Input data:
            - List of SampleIDs: (default location)`SampleIDs/igsr_samples.tsv.tsv`
            - Specify population IDs: `ASW` `CEU` `YRI`
        - Output (in SampleIDs folder):
            - `ASW_CEU_YRI_IDs.txt`
            - `ASW_Sample_IDs_haploid.txt`
            - `CEU_Sample_IDs_haploid.txt`
            - `YRI_Sample_IDs_haploid.txt`

                This file is all the sample IDs across the populations specified at the input.
                e.g. HG00318   HG00319    HH00320
                Each of these files has the sample IDs across a specific population, written as a haploid ( _1 and _2):
                 e.g. HG00318_1    HG00318_2     HG00319_1     HG00319_2

#### 2018 TODOs ####
- Right now the transition matrix is based on population-level ancestry proportions, so the transition matrix is the same for each individual. It would probably be more accurate if the transition matrix is calculated for each individual based on that individual's estimated ancestry proportions (estimated using ADMIXTURE). Currently, STRUCTUREpainter only expects the results of ADMIXTURE on the training set, not on the evaluation set ("canvas" set?), which means that this is not a trivial change.

#### Data Sources ####
- The file popgen/pops_for_sample_IDs.tsv can be downloaded from http://www.internationalgenome.org/data-portal/sample . Select the populations of interest (YRI, CEU and ASW) from the checkboxes, and then select the option to "Download the list" at the top of the page.
	- Rename the file as pops_for_sample_IDs.tsv
- The input data files for the vcftools command can be downloaded from http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ or ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
    - As mentioned on this page, the VCF files produced by the final phase of the 1000 Genomes Project (phase 3) are phased (http://www.internationalgenome.org/faq/are-1000-genomes-variant-calls-phased/)


3. Select out populations of interest: ASW, CEU, and YRI.

    In this script, we are taking only the phase3 data for ALL the populations, and extracting the phase3 data
     for the populations for which we extracted the IDs in the previous step, in (pops)_ID.txt file.
     The output file is a combined file across ALL the populations of interest. The breakup of this file into
     the individual populations is done in a later step.
     
    - Use the bash script: `bin/bashScripts/new/select_populations.sh`.
       - Script runtime for Chr22: about 8 minutes.
        - Calls:
            - vcftools
    - Input data:
        - In vcf folder: `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`
        - In SampleIDs: `ASW_CEU_YRI_IDs.txt`
    - Output: `chr22.phase3.ASW_CEU_YRI.SNPs.recode.vcf`, `chr22.phase3.ASW_CEU_YRI.SNPs.log`


3. Split the phased chromosomes into separate chromosomes

    This script just operates on the "combined" vcf file generated above
    - Script runtime for Chr22: about 10 minutes
    - Output: 
        - `chr22.ASW_CEU_YRI.SNPs.homologous.25k.2.Q`
        - `chr22.ASW_CEU_YRI.SNPs.homologous.25k.2.P`
        - `chr22.phase3.ASW_CEU_YRI.SNPs.homologous.vcf`
        - `chr22.ASW_CEU_YRI.SNPs.homologous.header.vcf`
        - `chr22.ASW_CEU_YRI.SNPs.homologous.25k.bed`
        - `chr22.phase3.ASW_CEU_YRI.SNPs.allele.vcf`
        # @Todo This should be allele.vcf.txt

#### Additional Info ####
As part of initial verification, admixture was tested using a third, real population which was the result of admixture between two other real populations. See docs for more details.

- New ToDos
    - @ToDo The vcftools command below is simply comverting the vcf file to a tped file format. Why do we need this if plink --vcf can accept a vcf file directly as input.
    - @ToDo Probably this is done because --vcf option is only available in plink 1.9 (not plink 1.07).
    - @ToDo Switch using plink 1.9 and eliminate the need for tped file generation
    - @ToDo Delete the file data
    - @ToDo Add documentation about verifying Admixture
    - @ToDo Figure out which Plink is being used and works
    - @ToDo Admixture and Admixed folders need to be moved within the {pops}_Data/Chr{num}/ folder. This way, each population set for each chromosome gets their own results for Structurepainter# STRUCTUREpainter
