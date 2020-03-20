# CRISPR_screen_analysis

A collection of python scripts to analyze data from single knockout or double knockout CRISPR genetic screens. 
I updated the readme file with brief descriptions how to execute the scripts for single knockout CRISPR screens (March 19, 2020)
I am planning to update this manual with more detailed explanations

Email Kyuho Han (kyuho@medic-life-sciences.com) for more info.

Copyright (c) 2020, Kyuho Han

Redistribution and use in source are permitted as long as the copyright notice above is included

If you use these scripts please cite the following paper :
Kyuho Han, et al., CRISPR screens in cancer spheroids identify 3D growth-specific vulnerabilities, Nature, 2020


###############################################################################
# Quick use:
###############################################################################

(1) Combine four fastq.gz files with the same header, but four different flow cell lane numbers
into single fastq.gz file. This works for all fastq.gz files in a designated folder

python GetSingleReadFastq.py <folder with demultiplexed fastq.gz files from bcl2fastq demultiplexing>

ex) python GetSingleReadFastq.py ../fastq/Batch_retest/


(2) Align fastq.gz files against an index file of a sgRNA library to generate count files

python GetSingleCounts.py <folder with fastq.gz files processed in step (1)> <output file> <header of bowtie index file in "indices" folder>

ex) python GetSingleCounts.py ../fastq/Batch_retest/ ../counts/Batch_retest/ Lung3D_Retest


(3) Compare two count files to calculate enrichment (log fold enrichment) scores of sgRNAs가

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