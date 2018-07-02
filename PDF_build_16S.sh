#!/usr/bin/sh
#######################################################################
mkdir 09PDFbuild
cd 09PDFbuild
mkdir figures
mkdir figures/0_raw_reads
mkdir figures/1_clean_reads
mkdir figures/0_raw_reads/0_1_base_quality
for file in `ls ../00_raw_reads/fastqc_out/*_fastqc/Images/per_base_quality.png`;do files=${file/..\/00_raw_reads\/fastqc_out\//};files=${files/_fastqc\/Images\/per_base_quality.png/};echo $files;cp $file figures/0_raw_reads/0_1_base_quality/$files.png;done
#######################################################################
mkdir figures/1_clean_reads/1_1_base_quality
for file in `ls ../01_clean_reads/fastqc_out/*_fastqc/Images/per_base_quality.png`;do files=${file/..\/01_clean_reads\/fastqc_out\//};files=${files/_fastqc\/Images\/per_base_quality.png/};echo $files;cp $file figures/1_clean_reads/1_1_base_quality/$files.png;done
#######################################################################
mkdir figures/2_OTU_assembly_analysis
mkdir figures/2_OTU_assembly_analysis/2_1_assemble_table 
mkdir figures/2_OTU_assembly_analysis/2_2_assemble_table
cp ../05_OTU_table_build/otu_tax_table.biom figures/2_OTU_assembly_analysis/
cp ../05_OTU_table_build/otu_tax_table_norm.biom figures/2_OTU_assembly_analysis/
cp ../05_OTU_table_build/otu_tax_table_norm.txt figures/2_OTU_assembly_analysis/otu_tax_table_norm.xls
cp ../05_OTU_table_build/otu_tax_table.txt figures/2_OTU_assembly_analysis/otu_tax_table.xls
cp ../03_pick_OTU/otus_table.txt figures/2_OTU_assembly_analysis/otus_table.xls
cp ../03_pick_OTU/rep_seq.fa figures/2_OTU_assembly_analysis/
cp ../03_pick_OTU/rep_seq_phylo.tre figures/2_OTU_assembly_analysis/
cp ../03_pick_OTU/Statistic_RawCleanNums.xls figures/2_OTU_assembly_analysis/2_2_assemble_table/sample2sequence_summary.xls
cp ../05_OTU_table_build/*.xls figures/2_OTU_assembly_analysis/
#########################################################################
mkdir figures/3_rarefaction
cp ../06_rare_analysis/*.svg figures/3_rarefaction/
##################################################
