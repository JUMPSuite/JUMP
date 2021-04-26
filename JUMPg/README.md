----------------------
Introduction
----------------------

JUMPg is a proteogenomics software pipeline for analyzing large mass spectrometry (MS) and functional genomics datasets. The pipeline includes customized database building, tag-based database search, peptide-spectrum match filtering, and data visualization. The program is written in perl and is designed for high performance parallel computing systems. This version of JUMPg depends on the proteomics software JUMP suite (https://github.com/JUMPSuite/JUMP).
 
For more information, please read: Li et al., JUMPg: an integrative proteogenomics pipeline identifying unannotated proteins in human brain and cancer cells. J. Proteome Res., 15 (2016), pp. 2309–2320

----------------------
Contents of this file
----------------------

* Introduction
* Software Requirements
* Hardware Requirements
* Installation
* Command Line Arguments
* Maintainers

----------------------
Software Requirements
---------------------- 

This version of JUMPg depends on the installation of proteomics software JUMP suite (https://github.com/JUMPSuite/JUMP).

----------------------
Hardware Requirements
---------------------- 

The program has been successfully tested on the following system:

+ Cluster mode (key parameters: 'cluster = 1' & 'Job_Management_System = SGE'):
  - Batch-queuing system: SGE, version 6.1u5, qsub and qstat
  - 128 GB memory and 64 cores on each node

+ Single server mode (key parameters: 'cluster = 0' & 'processors_used = 8'): 
  - 32 GB memory
  - 2 GHz CPU processors with 8 cores
 
----------------------
Installation
---------------------- 

INSTALLATION and CONFIGURATION:

__step 0__: install JUMP suite: 
 - follow instructions on JUMP suite GitHub site: https://github.com/JUMPSuite/JUMP

__Step 1__: download AnnoVar annotation files by running the following commands (take human hg19 as an example):
 - mkdir annotations
 - cd annotations
 - perl ../programs/c/customizedDB/annovar/annotate_variation.pl -downdb -buildver hg19 -downdb knownGene human/
 - cd ..
 
__Step 2__: prepare reference genome (FASTA format)  # this file is required for splice junction and RNA 6FT analysis
 - For model organisms, the reference genome should be available through UCSC genome browser database (http://genome.ucsc.edu/). 
 
   Example command: 
   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
   Concatenate all the fasta files using the following command
   cat *.fa > hg19_genome_raw.fa

----------------------
Run the example
----------------------

Since the program supports multistage database search analysis, we have designed a two-stage test dataset. For the 1st stage, the MS/MS spectra are searched against a peptide database pooled from uniProt, mutations and splice junctions; the matching results are filtered to ~1% false discovery rate. For the 2nd stage, the remaining high quality spectra are searched against the 6FT database.

Users may download this example dataset and the asscoiated resources (e.g., databases and human genome annotation files) from the JUMPg GitHub site: https://github.com/liyuxin-bioinformatics/JUMPg/tree/master/exampleData_v6.2

__1st stage analysis:__

__Step 1__: cp 'jump_g_v2.2.stage1.params' (included in exampleData_v6.2/parameters folder) to your working directory and edit the following parameters:

The following assumes that your working directory is at "/home/usr/JUMPg"

 - input_ref_proteins = /home/usr/JUMPg/exampleData_v6.2/uniProt/reference_uniProt_human.fasta
 - input_mutation = /home/usr/JUMPg/exampleData_v6.2/rna_database/mutations.txt
 - input_junction = /home/usr/JUMPg/exampleData_v6.2/rna_database/junctions.txt
 - annovar_annotation_file = /home/usr/JUMPg/JUMPg_v2.3.1/annotations/human/hg19_knownGene.txt
 - gene_ID_convert_file = /home/usr/JUMPg/JUMPg_v2.3.1/annotations/human/hg19_kgXref.txt
 - reference_genome = /home/usr/genomes/hg19_genome_raw.fa

__Step 2__: cp spectra file 'exampleData_v6.2/spectra/MS_test.mzXML' (included in exampleData_v6.2.tar.gz) to your working diredirectory and run the program using the following command:

> perl /home/usr/JUMPg_v2.3.1/programs/JUMPg_v2.3.1.pl jump_g_v2.2.stage1.params MS_test.mzXML

__Output__: the program will generate a folder named 'gnm_round1_test1'. The final results are all included in the 'publications' folder that includes 6 files:
 - identified_peptide_table.txt: a text table, showing the identified peptide sequences, PSM counts, tags and scores of the best PSM, and database entries.
 - mutation_peptides.txt: a text table showing the identified mutation peptides
 - mutation_peptides.bed: mutation peptides with genomic locations in BED format, which can be co-displayed with other genomic information on the UCSC genome browser
 - junction_peptides.txt: a text table showing the identified novel junction peptides
 - junction_peptides.bed: novel junction peptides with genomic locations in BED format for visualization on the UCSC genome browser
 - reference_peptides.bed: reference peptides with genomic locations in BED format for visualization on the UCSC genome browser

For multistage analysis, the program also generates the unmatched high quality MS/MS spectra, of which the path is recorded in this file: gnm_round1_test1/multistage/qc_MSMS_input.txt. This file will be taken as input for 2nd stage analysis.

__2nd stage analysis:__

__Step 1__: cp 'jump_g_v2.2.stage2.params' (included in exampleData_v6.2/parameters) to your working directory and edit the following parameters:
 - input_transcript_6FT = /home/usr/JUMPg/exampleData_v6.2/rna_database/assembled_transcripts.fasta
 - annovar_annotation_file = /home/usr/JUMPg/JUMPg_v2.3.1/annotations/human/hg19_knownGene.txt
 - gene_ID_convert_file = /home/usr/JUMPg/JUMPg_v2.3.1/annotations/human/hg19_kgXref.txt
 - reference_genome = /home/usr/genomes/hg19_genome_raw.fa
 - reference_protein = /home/usr/JUMPg/exampleData_v6.2/uniProt/reference_uniProt_human.fasta

__Step 2__: copy qc_MSMS_input.txt (that records the path to unmatched high quality MS/MS spectra) to current directory:

>cp gnm_stage1_test1/multistage/qc_MSMS_input.txt .

__Step 3__: run the program by the command:
>perl /home/usr/JUMPg_v2.3.1/programs/JUMPg_v2.3.1.pl jump_g_v2.2.stage2.params qc_MSMS_input.txt

__Output__: similar to the 1st stage result, the program will generate a folder named 'gnm_round2_test1' containing results in its 'publications' folder.


----------------------
Maintainers
----------------------

* To submit bug reports and feature suggestions, please contact:

  Yuxin Li (yuxin.li@stjude.org) and Junmin Peng (junmin.peng@stjude.org)

