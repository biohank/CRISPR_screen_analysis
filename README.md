# CRISPR_screen_analysis

A collection of python scripts to analyze data from single knockout or double knockout CRISPR genetic screens. 
I updated the readme file with brief descriptions how to execute the scripts for single knockout CRISPR screens (March 19, 2020)
I am planning to update this manual with more detailed explanations

Copyright (c) 2020, Kyuho Han

Redistribution and use in source are permitted as long as the copyright notice above is included

If you use these scripts please cite the following paper :
Kyuho Han, et al., CRISPR screens in cancer spheroids identify 3D growth-specific vulnerabilities, Nature, 2020


###############################################################################
# Quick use for single sgRNA CRISPR library
###############################################################################

(1) Combine four fastq.gz files with the same header into a single file in a designated folder

	python GetSingleReadFastq.py <folder with demultiplexed fastq.gz files> <output folder for combined fastq>

	ex) python GetSingleReadFastq.py ../fastq/Batch_retest/ ../fastq/Batch_retest_all/

(2) Align fastq.gz files against an index file of a sgRNA library to generate count files

    python GetSingleCounts.py <folder with fastq.gz files processed in step (1)> <output file> <header of bowtie index file in "indices" folder>

    ex) python GetSingleCounts.py ../fastq/Batch_retest/ ../counts/Batch_retest/ Lung3D_Retest


(3) Compare two count files to calculate enrichment (log fold enrichment) scores of sgRNAs

	python GetEnrichment.py <count file 1> <count file 2> <outputfile> -gt single -ct <threshold for sgRNA count>

	ex) python GetEnrichment.py ../counts/Batch_retest/T0_23_r1.counts ../counts/Batch_retest/Day21_23_rep1_r1.counts ../results/Batch_retest/T0_vs_Day21_23_rep1 -gt single -ct 10
	ex) python GetEnrichment.py ../counts/Batch_retest/T0_23_r1.counts ../counts/Batch_retest/Day21_23_rep2_r1.counts ../results/Batch_retest/T0_vs_Day21_23_rep2 -gt single -ct 10


(4) Calculate gene phenotypes (calculate the median enrichment score of all sgRNAs that target a gene)

	python GetPhenotype.py <sgRNA enrichment csv file>

	ex) python GetPhenotype.py ../results/Batch_retest/T0_vs_Day21_23_rep1
	ex) python GetPhenotype.py ../results/Batch_retest/T0_vs_Day21_23_rep2


(5) Combine two replicates to calculate gene phenotypes

	python CombineReplicates.py <sgRNA enrichment csv file> <output file>

	ex) python CombineReplicates.py ../results/Batch_retest/T0_vs_Day21_23_rep1 ../results/Batch_retest/T0_vs_Day21_23_rep2 ../results/Batch_retest/T0_vs_Day21_23_combo


###############################################################################
# Quick use for double sgRNA CRISPR library
###############################################################################

(1) Generate count files for double sgRNAs using GetDoubleCounts.py

python GetDoubleCounts.py <folder that contains all fastq.gz files> <output folder for count files> <index name for front sgRNAs> <index name for rear sgRNAs>

ex1) python GetDoubleCounts.py /mnt/lab_data/bassik/kyuhohan/NextSeq/bcl2fastq/KRAS_PPI_20by20_180206_Y_Y_I/ /mnt/lab_data/bassik/kyuhohan/Kyuho_Screening_Analysis/counts/20180214_KRAS_PPI_20by20/ KRAS_20by20 KRAS_20by20
ex2) python GetDoubleCounts.py /mnt/lab_data/bassik/kyuhohan/NextSeq/bcl2fastq/KRAS_pgRNA_Kaitlyn_180316/ /mnt/lab_data/bassik/kyuhohan/Kyuho_Screening_Analysis/counts/20180316_KRAS_PPI_pgRNA/ KRAS_PPI_pgRNA KRAS_PPI_pgRNA 

This will generate count files which you can use in the next script

(2) Compare two count files to calculate enrichment (log fold enrichment) scores of double sgRNAs using GetEnrichment.py

python GetEnrichment.py <count file 1> <count file 2> <outputfile header>

ex1) python GetEnrichment.py ../counts/20171009_SL_Screens/20171009_SL_Plas_pgDouble_Final.counts ../counts/20171009_SL_Screens/20171009_SL_D15_Unt1_pgDouble_Final.counts ../results/20171009_SL_Screens/Plas_vs_D15_Unt1

ex2) python GetEnrichment.py ../counts/20171009_SL_Screens/20171009_SL_Plas_pgDouble_Final.counts ../counts/20171009_SL_Screens/20171009_SL_D15_Unt2_pgDouble_Final.counts ../results/20171009_SL_Screens/Plas_vs_D15_Unt2

This will generate enrichment.csv files which you can use in the next script

(3) Calculate gene phenotypes (calculate the median enrichment score of all double sgRNAs that target a gene pair) using GetPhenotype.py

python GetPhenotype.py <header of sgRNA enrichment csv file>

ex1) python GetPhenotype.py ../results/20171009_SL_Screens/Plas_vs_D15_Unt1
ex2) python GetPhenotype.py ../results/20171009_SL_Screens/Plas_vs_D15_Unt2

This will generate phenotype.csv files which you can use in the next script

(4) Combine two replicates to calculate gene phenotypes using CombineReplicates.py

python CombineReplicates.py <phenotype.csv file from replicate 1> <phenotype.csv file from replicate 2> <outputfile header for combo file>

ex1) python CombineReplicates.py ../results/20171009_SL_Screens/Plas_vs_D15_Unt1 ../results/20171009_SL_Screens/Plas_vs_D15_Unt2 ../results/20171009_SL_Screens/Plas_vs_D15_Unt_Combo

This will generate final combo files
