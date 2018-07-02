#!/usr/bin/sh
#####use it in 07_group/group_*****/
cd "07_groups/group_"$1
mkdir need_file_folder need_file_folder/4_alpha_diversity need_file_folder/5_Anosim_analysis need_file_folder/6_beta_diversity need_file_folder/7_environment_factor need_file_folder/8_species_visulization need_file_folder/9_differential_species_analysis 
mkdir need_file_folder/10_function_analysis
mkdir need_file_folder/6_beta_diversity/6_1_ordination need_file_folder/6_beta_diversity/6_2_heatmap need_file_folder/6_beta_diversity/6_3_clust
mkdir need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_without_clust
mkdir need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_with_clust
mkdir need_file_folder/6_beta_diversity/6_3_clust/6_3_1_no_bar
mkdir need_file_folder/9_differential_species_analysis/9_1_LEfSe_analysis 
mkdir need_file_folder/9_differential_species_analysis/9_1_LEfSe_analysis/differential_species_distribution

#############################################################################
cp 09_alpha_analysis/*.svg need_file_folder/4_alpha_diversity/richness_plot.svg
cp 09_alpha_analysis/*.txt need_file_folder/4_alpha_diversity/multi_adiversity.xls
##############################################################################
cp 08_Anosim_analysis/*.svg  need_file_folder/5_Anosim_analysis/anosim_plot.svg
cp otu_tab_norm_* need_file_folder/5_Anosim_analysis/
rm need_file_folder/5_Anosim_analysis/*.biom
mv need_file_folder/5_Anosim_analysis/*.txt need_file_folder/5_Anosim_analysis/otu_table_norm.xls
##############################################################################
cp 10_beta_analysis/pca* need_file_folder/6_beta_diversity/6_1_ordination/pca_plot.svg
cp 10_beta_analysis/NMDS* need_file_folder/6_beta_diversity/6_1_ordination/NMDS_plot.svg
cp 10_beta_analysis/PCoA* need_file_folder/6_beta_diversity/6_1_ordination/PCoA_plot.svg
####################
cp 11_tax_classification/otu_taxlevel_2.txt need_file_folder/6_beta_diversity/6_2_heatmap/phylum_tax.xls
cp 11_tax_classification/otu_taxlevel_3.txt need_file_folder/6_beta_diversity/6_2_heatmap/class_tax.xls
cp 11_tax_classification/otu_taxlevel_4.txt need_file_folder/6_beta_diversity/6_2_heatmap/order_tax.xls
cp 11_tax_classification/otu_taxlevel_5.txt need_file_folder/6_beta_diversity/6_2_heatmap/family_tax.xls
cp 11_tax_classification/otu_taxlevel_6.txt need_file_folder/6_beta_diversity/6_2_heatmap/genus_tax.xls
cp 11_tax_classification/otu_taxlevel_7.txt need_file_folder/6_beta_diversity/6_2_heatmap/species_tax.xls
#####
#cp 11_tax_classification/heatmap* need_file_folder/6_beta_diversity/6_2_heatmap/mv need_file_folder/6_beta_diversity/6_2_heatmap/*noclust.svg need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_without_clust/
#mv need_file_folder/6_beta_diversity/6_2_heatmap/*.svg need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_with_clust/
cp 11_tax_classification/heatmap_level_2_noclust.svg need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_without_clust/heatmap_p_noclust.svg
cp 11_tax_classification/heatmap_level_3_noclust.svg need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_without_clust/heatmap_c_noclust.svg
cp 11_tax_classification/heatmap_level_4_noclust.svg need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_without_clust/heatmap_o_noclust.svg
cp 11_tax_classification/heatmap_level_5_noclust.svg need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_without_clust/heatmap_f_noclust.svg
cp 11_tax_classification/heatmap_level_6_noclust.svg need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_without_clust/heatmap_g_noclust.svg
cp 11_tax_classification/heatmap_level_7_noclust.svg need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_without_clust/heatmap_s_noclust.svg
###
 cp 11_tax_classification/heatmap_level_2.svg need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_with_clust/heatmap_p.svg
 cp 11_tax_classification/heatmap_level_3.svg need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_with_clust/heatmap_c.svg
 cp 11_tax_classification/heatmap_level_4.svg need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_with_clust/heatmap_o.svg
 cp 11_tax_classification/heatmap_level_5.svg need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_with_clust/heatmap_f.svg
 cp 11_tax_classification/heatmap_level_6.svg need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_with_clust/heatmap_g.svg
 cp 11_tax_classification/heatmap_level_7.svg need_file_folder/6_beta_diversity/6_2_heatmap/6_2_1_with_clust/heatmap_s.svg

#####
cp 11_tax_classification/otu_taxlevel_2_clust.svg need_file_folder/6_beta_diversity/6_3_clust/6_3_1_no_bar/phylum_cluster.svg
cp 11_tax_classification/otu_taxlevel_3_clust.svg need_file_folder/6_beta_diversity/6_3_clust/6_3_1_no_bar/class_cluster.svg
cp 11_tax_classification/otu_taxlevel_4_clust.svg need_file_folder/6_beta_diversity/6_3_clust/6_3_1_no_bar/order_cluster.svg
cp 11_tax_classification/otu_taxlevel_5_clust.svg need_file_folder/6_beta_diversity/6_3_clust/6_3_1_no_bar/family_cluster.svg
cp 11_tax_classification/otu_taxlevel_6_clust.svg need_file_folder/6_beta_diversity/6_3_clust/6_3_1_no_bar/genus_cluster.svg
cp 11_tax_classification/otu_taxlevel_7_clust.svg need_file_folder/6_beta_diversity/6_3_clust/6_3_1_no_bar/species_cluster.svg

cp 11_tax_classification/otu_taxlevel_2_clust_bar.svg need_file_folder/6_beta_diversity/6_3_clust/phylum_clust_bar.svg
cp 11_tax_classification/otu_taxlevel_3_clust_bar.svg need_file_folder/6_beta_diversity/6_3_clust/class_clust_bar.svg
cp 11_tax_classification/otu_taxlevel_4_clust_bar.svg need_file_folder/6_beta_diversity/6_3_clust/order_clust_bar.svg
cp 11_tax_classification/otu_taxlevel_5_clust_bar.svg need_file_folder/6_beta_diversity/6_3_clust/family_clust_bar.svg
cp 11_tax_classification/otu_taxlevel_6_clust_bar.svg need_file_folder/6_beta_diversity/6_3_clust/genus_clust_bar.svg
cp 11_tax_classification/otu_taxlevel_7_clust_bar.svg need_file_folder/6_beta_diversity/6_3_clust/species_clust_bar.svg
###############################3
cp need_file_folder/6_beta_diversity/6_2_heatmap/*.xls need_file_folder/6_beta_diversity/6_3_clust/
########
cp 14_species_tree_distribution/*.svg need_file_folder/8_species_visulization/
cp otu_tab_norm* need_file_folder/8_species_visulization/
cp 11_tax_classification/otu_taxlevel_2.svg need_file_folder/8_species_visulization/phylum_bar.svg
cp 11_tax_classification/otu_taxlevel_3.svg need_file_folder/8_species_visulization/class_bar.svg
cp 11_tax_classification/otu_taxlevel_4.svg need_file_folder/8_species_visulization/order_bar.svg
cp 11_tax_classification/otu_taxlevel_5.svg need_file_folder/8_species_visulization/family_bar.svg
cp 11_tax_classification/otu_taxlevel_6.svg need_file_folder/8_species_visulization/genus_bar.svg
cp 11_tax_classification/otu_taxlevel_7.svg need_file_folder/8_species_visulization/species_bar.svg
##############
cp 12_LEfSe_analysis/*.svg need_file_folder/9_differential_species_analysis/9_1_LEfSe_analysis
cp 12_LEfSe_analysis/*.xls need_file_folder/9_differential_species_analysis/9_1_LEfSe_analysis
cp 12_LEfSe_analysis/biomarkers*/* need_file_folder/9_differential_species_analysis/9_1_LEfSe_analysis/differential_species_distribution/
##############
#cp 11_tax_classification/otu_taxlevel_2.svg ..//..//need_file/2_OTU_assembly_analysis/2_1_assemble_table/phylum_barplot.svg
#cp 11_tax_classification/otu_taxlevel_3.svg ..//..//need_file/2_OTU_assembly_analysis/2_1_assemble_table/class_barplot.svg
#cp 11_tax_classification/otu_taxlevel_4.svg ..//..//need_file/2_OTU_assembly_analysis/2_1_assemble_table/order_barplot.svg
#cp 11_tax_classification/otu_taxlevel_5.svg ..//..//need_file/2_OTU_assembly_analysis/2_1_assemble_table/family_barplot.svg
#cp 11_tax_classification/otu_taxlevel_6.svg ..//..//need_file/2_OTU_assembly_analysis/2_1_assemble_table/genus_barplot.svg
#cp 11_tax_classification/otu_taxlevel_7.svg ..//..//need_file/2_OTU_assembly_analysis/2_1_assemble_table/species_barplot.svg
#cp 11_tax_classification/otu_taxlevel_2.txt ..//..//need_file/2_OTU_assembly_analysis/2_1_assemble_table/phylum_tax.xls
#cp 11_tax_classification/otu_taxlevel_3.txt ..//..//need_file/2_OTU_assembly_analysis/2_1_assemble_table/class_tax.xls
#cp 11_tax_classification/otu_taxlevel_4.txt ..//..//need_file/2_OTU_assembly_analysis/2_1_assemble_table/order_tax.xls
#cp 11_tax_classification/otu_taxlevel_5.txt ..//..//need_file/2_OTU_assembly_analysis/2_1_assemble_table/family_tax.xls
#cp 11_tax_classification/otu_taxlevel_6.txt ..//..//need_file/2_OTU_assembly_analysis/2_1_assemble_table/genus_tax.xls
#cp 11_tax_classification/otu_taxlevel_7.txt ..//..//need_file/2_OTU_assembly_analysis/2_1_assemble_table/species_tax.xls
mkdir -p need_file_folder/10_function_analysis/10_1_KEGG
mkdir -p need_file_folder/10_function_analysis/10_2_COG
cp "../../08_Function_predict/group_"$1"/kegg_bar_"$1".svg" need_file_folder/10_function_analysis/10_1_KEGG/KEGG_pathway.svg
cp "../../08_Function_predict/group_"$1"/ko_tab_input_"$1".txt" need_file_folder/10_function_analysis/10_1_KEGG/kegg_pathway_prediction_count.xls
cp "../../08_Function_predict/group_"$1"/cog_tab_input_"$1".txt" need_file_folder/10_function_analysis/10_2_COG/cog_prediction_count.xls
cp "../../08_Function_predict/group_"$1"/cog_bar_"$1".svg" need_file_folder/10_function_analysis/10_2_COG/COG_function.svg
cp "../../08_Function_predict/group_"$1"/ANOVAtab_cog_"$1".xls" need_file_folder/10_function_analysis/10_2_COG/
cp "../../08_Function_predict/group_"$1"/ANOVAtab_kegg_"$1".xls" need_file_folder/10_function_analysis/10_1_KEGG/
cp ../../08_Function_predict/picrust_otus_css.biom  need_file_folder/10_function_analysis/
mkdir -p need_file_folder/theme
inkscape -f 08_Anosim_analysis/Anosim_plot*.svg  -z -d 300 -e need_file_folder/theme/ANOSIM.png
inkscape -f 09_alpha_analysis/*.svg -z -d 300 -e need_file_folder/theme/alpha_diversity_boxplot.png
inkscape -f 10_beta_analysis/NMDS*.svg -z -d 300 -e need_file_folder/theme/nmds分析.png
inkscape -f 10_beta_analysis/pca_plot*.svg  -z -d 300 -e need_file_folder/theme/pca分析.png
inkscape -f 10_beta_analysis/PCoA*.svg -z -d 300 -e need_file_folder/theme/pcoa分析.png
inkscape -f 11_tax_classification/otu_taxlevel_6.svg  -z -d 300 -e need_file_folder/theme/样本菌群结构.png
inkscape -f 11_tax_classification/otu_taxlevel_6_clust_bar.svg -z -d 300 -e need_file_folder/theme/样本菌群结构聚类.png
inkscape -f 11_tax_classification/heatmap_level_6.svg -z -d 300 -e need_file_folder/theme/heatmap_f.png
inkscape -f 11_tax_classification/heatmap_level_6_noclust.svg -z -d 300 -e need_file_folder/theme/heatmap_f_noclust.png
inkscape -f 11_tax_classification/otu_taxlevel_6_clust.svg -z -d 300 -e need_file_folder/theme/genus_cluster.png
inkscape -f 12_LEfSe_analysis/LEfSe_caladogram_*.svg -z -d 300 -e need_file_folder/theme/LEfSe分析2.png
inkscape -f 12_LEfSe_analysis/LEfSe_LDA*.svg -z -d 300 -e need_file_folder/theme/LEfSe分析.png
cp "12_LEfSe_analysis/biomarkers_"$1"/LEfSe_DET_Figure1.png" need_file_folder/theme/LEfSe分析3.png
for file in `ls ../../01_clean_reads/fastqc_out/*_clean_merge_fastqc/Images/per_base_quality.png`;do cp $file need_file_folder/theme/per_base_quality.png;done
inkscape -f Venn_plot*.svg -z -d 300 -e need_file_folder/theme/样品OTU分布Venn图.png
inkscape -f 14_species_tree_distribution/Species_Tree*.svg -z -d 300 -e need_file_folder/theme/样本菌群总览图.png
#inkscape -f ..//..//08_Function_predict/group_*/cog_metagenome_predict.OTA.svg -z -d 300 -e theme/cog.png
#inkscape -f ..//..//08_Function_predict/group_group/ko_metagenome_predict_L3.OTA.svg -z -d 300 -e theme/kegg.png
inkscape -f ../../06_rare_analysis/observed_species.svg -z -d 300 -e need_file_folder/theme/RarefactionCurve.png
inkscape -f ../../06_rare_analysis/shannon.svg -z -d 300 -e need_file_folder/theme/Shannon-Wiener曲线.png
inkscape -f ../../06_rare_analysis/simpson.svg -z -d 300 -e need_file_folder/theme/Simpson稀释曲线.png
inkscape -f "../../08_Function_predict/group_"$1"/kegg_bar_"$1".svg" -z -d 300 -e need_file_folder/theme/kegg.png
inkscape -f "../../08_Function_predict/group_"$1"/cog_bar_"$1".svg" -z -d 300 -e need_file_folder/theme/cog.png

mv need_file_folder/ "../../09PDFbuild/figures/group_"$1
