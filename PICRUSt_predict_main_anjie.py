#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import argparse
def readparas(args):
	parser = argparse.ArgumentParser(description="The pipeline was designed to run the ITS and 16s analysis and visulization.\n"
		"The developer: ShuangbinXu.\n"
		"Email:xusbin@anjiemed.com\n")
	parser.add_argument('configure_file',metavar='configure_file', type=str, help="the configure file.")
	args = parser.parse_args()
	return vars(args)

def congfile(configurefile):
	filedict = {}
	with open(configurefile, 'rt') as fin:
		for line in fin:
			tmp = line.replace("\n", "").split("=")
			filedict[tmp[0]] = tmp[1]
	return filedict

def mkdir_PICRUSt_dir(workpath):
	#if not os.path.exists(workpath+"predict_function_dir.sh"):
	with open("./predict_function_dir.sh", "w") as dirfile:
		if not os.path.exists(workpath+"/Shellscript"):
			dirfile.write("#!/bin/bash\ncd " + workpath + "\nmkdir 08_Function_predict Shellscript")
			dirfile.close()
		else:
			dirfile.write("#!/bin/bash\nmkdir 08_Function_predict")
			dirfile.close()
		if not os.path.exists(workpath+"/08_Function_predict"):
			os.system("sh predict_function_dir.sh")

def subgroup_PICRUSTdir(workpath, groups_list, flaggroup):
	tmpgroup = groups_list.split(";")
	tmp = []
	grouplist = []
	for i in xrange(len(tmpgroup)):
		groupes = tmpgroup[i].split(",")
		tmpdir = workpath+"/08_Function_predict/group_"+"_vs_".join(groupes)
		tmp.append("mkdir -p " + tmpdir )#+ "\ncd "+tmpdir)
	groupsh = "\n".join(tmp)
	if flaggroup=="group":
		with open(workpath+"/Shellscript/subgroupdir1_picrust.sh", "w") as subdirfile:
			subdirfile.write("#!/usr/bin/sh\n"+groupsh)
			subdirfile.close()
			if not os.path.exists(workpath+"/08_Function_predict/group_"+"_vs_".join(groupes)) and os.path.exists(workpath+"/07_groups/"):
				os.system("sh "+workpath+"/Shellscript/subgroupdir1_picrust.sh")
	if flaggroup=="subgroup":
		with open(workpath+"/Shellscript/subgroupdir2_picrust.sh", "w") as subdirfile:
			subdirfile.write("#!/usr/bin/sh\n"+groupsh)
			subdirfile.close()
			if not os.path.exists(workpath+"/08_Function_predict/group_"+"_vs_".join(groupes)) and os.path.exists(workpath+"/07_groups/"):
				os.system("sh "+workpath+"/Shellscript/subgroupdir2_picrust.sh")

def picrust_otu(workpath, tax_fasta, taxonomy, runPickOTU):
	if not os.path.exists(workpath+"/08_Function_predict/08_pick_close_OTU.log") and os.path.exists(workpath+"/02_assemble/usearch_input_merge.fa"):	
		with open(workpath+"/Shellscript/08_picrust_otu.sh" , "w") as shfile:
			shfile.write("#!/usr/bin/sh\ncd " + workpath + "/03_pick_OTU\nusearch10 -usearch_global ../02_assemble/usearch_input_merge.fa -db " + tax_fasta + " -otutabout otu_table_usearch_close.txt -strand both -id 0.97\nusearch_close_tax.py " + taxonomy +" otu_table_usearch_close.txt otu_tax_table_usearch_close.txt\nbiom convert -i otu_tax_table_usearch_close.txt -o otu_tax_table_usearch_close.biom --table-type=\"Taxon table\" --to-hdf5 --process-obs-metadata taxonomy\nnormalize_table.py -i otu_tax_table_usearch_close.biom -o "+workpath+"/08_Function_predict/picrust_otus_css.biom -a CSS\nnormalize_by_copy_number.py -i "+workpath+"/08_Function_predict/picrust_otus_css.biom -o " + workpath + "/08_Function_predict/normalized_picrust_otus.biom\npredict_metagenomes.py -i " + workpath + "/08_Function_predict/normalized_picrust_otus.biom -t cog -f -o " + workpath + "/08_Function_predict/cog_predict.txt\npredict_metagenomes.py -i " + workpath + "/08_Function_predict/normalized_picrust_otus.biom -t ko -o " + workpath + "/08_Function_predict/ko_predict.biom\ncategorize_by_function.py -i " + workpath + "/08_Function_predict/ko_predict.biom -c KEGG_Pathways -l 3 -f -o "+ workpath +"/08_Function_predict/ko_L3_predict.txt\ntouch "+workpath+"/08_Function_predict/08_pick_close_OTU.log")
			shfile.close()
			if runPickOTU == "Yes":
				os.system("sh "+workpath+"/Shellscript/08_picrust_otu.sh")

def splitPICRUStOTUtab(workpath, groups, subanalysis):
	tmpvs = "_vs_".join(groups.split(","))
	#subCOGTab = "cog_tab_"+tmpvs+".biom"
	subCOGTXT = "cog_tab_input_"+tmpvs+".txt"
	subCOGtxt = "cog_tab_"+tmpvs+".txt"
	subKOtxt = "ko_tab_"+tmpvs+".txt"
	subKOTXT = "ko_tab_input_"+tmpvs+".txt"
	flagCOG =  workpath+"/08_Function_predict/group_"+tmpvs+"/"+subCOGTXT
	flagKO = workpath+"/08_Function_predict/group_"+tmpvs+"/"+subKOtxt
	fileCOGtxt = workpath+"/08_Function_predict/group_"+tmpvs+"/"+subCOGtxt
	fileKOtxt = workpath+"/08_Function_predict/group_"+tmpvs+"/"+subKOtxt
	tmpsamp = workpath+"/Tmpgroups/group_"+tmpvs+"/subsample_"+tmpvs+".txt"
	fileCOGTXT = workpath+"/08_Function_predict/group_"+tmpvs+"/"+subCOGTXT
	fileKOTXT = workpath+"/08_Function_predict/group_"+tmpvs+"/"+ subKOTXT
	if os.path.exists(workpath+"/08_Function_predict/cog_predict.txt") and not os.path.exists(flagCOG):
		tmpsh = workpath+"/Tmpgroups/group_"+tmpvs+"/08_split_PICRUSt_"+tmpvs+".sh"
		with open(tmpsh, "w") as shfile:
			shfile.write("#!/usr/bin/sh\nextract_sample_otu_table.R "+workpath+"/08_Function_predict/cog_predict.txt "+ tmpsamp +" -o " + fileCOGtxt + "\nextract_sample_otu_table.R "+workpath+"/08_Function_predict/ko_L3_predict.txt "+tmpsamp+" -o " + fileKOtxt + "\nPICRUSt_out_change.py " + fileCOGtxt + " " + fileCOGTXT + " cog\nPICRUSt_out_change.py " + fileKOtxt +" " + fileKOTXT + " ko")
			shfile.close()
			if subanalysis == "Yes":
				os.system("sh "+tmpsh)

def PICRUSt_visualization(workpath, groups, subanalysis):
	tmpvs = "_vs_".join(groups.split(","))
	subCOGTXT = "cog_tab_input_"+tmpvs+".txt"
	subKOTXT = "ko_tab_input_"+tmpvs+".txt"
	fileCOGTXT = workpath+"/08_Function_predict/group_"+tmpvs+"/"+subCOGTXT
	fileKOTXT = workpath+"/08_Function_predict/group_"+tmpvs+"/"+ subKOTXT
	tmpsamp = workpath+"/Tmpgroups/group_"+tmpvs+"/subsample_"+tmpvs+".txt"
	koimage = workpath+"/08_Function_predict/group_"+tmpvs+"/kegg_bar_"+tmpvs+".svg"
	cogimage = workpath+"/08_Function_predict/group_"+tmpvs+"/cog_bar_"+tmpvs+".svg"
	annovatabcog = workpath+"/08_Function_predict/group_"+tmpvs+"/ANNOVAtab_cog_"+tmpvs+".svg"
	annovatabkegg = workpath+"/08_Function_predict/group_"+tmpvs+"/ANNOVAtab_kegg_"+tmpvs+".svg"
	if os.path.exists(fileCOGTXT) and not os.path.exists(koimage):
		tmpsh = workpath+"/Tmpgroups/group_"+tmpvs+"/08_PICRUSt_plot_"+tmpvs+".sh"
		with open(tmpsh, "w") as shfile:
			shfile.write("#!/usr/bin/sh\ncog_bar_plot.R "+fileCOGTXT+" "+tmpsamp+" "+cogimage+"\nkegg_path_bar_plot.R "+fileKOTXT+" "+ tmpsamp +" "+koimage+"\nANOVA_test_table.R "+fileCOGTXT+" "+tmpsamp+" -o "+annovatabcog+"\nANOVA_test_table.R "+fileKOTXT+" "+ tmpsamp +" -o "+annovatabkegg)			
			shfile.close()
			if subanalysis == "Yes":
				os.system("sh "+tmpsh)
	
			
if __name__ == "__main__":
	params = readparas(sys.argv)
	configure = params['configure_file']
	confiles = congfile(configure)
	workpath = confiles['workpath']
	PICRUSt_tax_fasta = confiles["tax_fasta"]
	PICRUSt_taxonomy = confiles['taxonomy']
	runPickOTU = confiles['runPickOTU']
	subanalysis = confiles["subanalysis"]
	groups = confiles["groups"]
	subgroups = confiles["subgroups"]
	mkdir_PICRUSt_dir(workpath)
	picrust_otu(workpath, PICRUSt_tax_fasta, PICRUSt_taxonomy, runPickOTU)
	if len(subgroups.split(";")[0]) != 0:
		subgroup_PICRUSTdir(workpath, subgroups, "subgroup")
	if len(groups.split(";")[0]) != 0:
		subgroup_PICRUSTdir(workpath, groups, "group")
	if len(groups.split(";")[0]) != 0:
		groupslist = groups.split(";")
		for i in xrange(len(groupslist)):
			splitPICRUStOTUtab(workpath, groupslist[i], subanalysis)
			PICRUSt_visualization(workpath, groupslist[i], subanalysis)
	if len(subgroups.split(";")[0]) != 0:
		subgroupslist = subgroups.split(";")
		for i in xrange(len(subgroupslist)):
			splitPICRUStOTUtab(workpath, subgroupslist[i], subanalysis)
			PICRUSt_visualization(workpath, subgroupslist[i], subanalysis)

