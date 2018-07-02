#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import os
import argparse
def readparas(args):
	parser = argparse.ArgumentParser(description="The pipeline was designed to run the ITS and 16s analysis and visulization.\n"
		"The developer: ShuangbinXu.\n"
		"Email:xusbin@anjiemed.com\n")
	parser.add_argument('sample_file',metavar="sample_file", type=str, help='the sample file')
	parser.add_argument('configure_file',metavar='configure_file', type=str, help="the configure file.")
	args = parser.parse_args()
	return vars(args)
def mkdir_16s():
	if not os.path.exists("dir_16.sh"):
		with open("./dir_16.sh", "w") as dirfile:
			dirfile.write("#!/bin/bash\nmkdir Shellscript 00_raw_reads 01_clean_reads 02_assemble 03_pick_OTU 04_OTU_assignment 05_OTU_table_build 06_rare_analysis 07_groups Tmpgroups\nmkdir -p 00_raw_reads/fastqc_out\nmkdir -p 01_clean_reads/fastqc_out")
			dirfile.close()
			os.system("sh dir_16.sh")

def subgroupdir(workpath,groups_list, flaggroup):
	tmpgroup = groups_list.split(";")
	tmp = []
	grouplist = []
	for i in xrange(len(tmpgroup)):
		groupes = tmpgroup[i].split(",")
		tmpdir = workpath+"/07_groups/group_"+"_vs_".join(groupes)
		tmpdirtmp = workpath+"/Tmpgroups/group_"+"_vs_".join(groupes)
		tmp.append("mkdir -p " + tmpdirtmp )
		tmp.append("mkdir -p " + tmpdir + "\ncd "+tmpdir)
		tmp.append("mkdir 08_Anosim_analysis 09_alpha_analysis 10_beta_analysis 11_tax_classification 12_LEfSe_analysis 13_environment_factor_analysis 14_species_tree_distribution")
	groupsh = "\n".join(tmp)
	if flaggroup=="group":
		with open(workpath+"/Shellscript/subgroupdir1.sh", "w") as subdirfile:
			subdirfile.write("#!/usr/bin/sh\n"+groupsh)
			subdirfile.close()
			if not os.path.exists(workpath+"/07_groups/group_"+"_vs_".join(groupes)) and os.path.exists(workpath+"/07_groups/"):
				os.system("sh "+workpath+"/Shellscript/subgroupdir1.sh")	
	if flaggroup=="subgroup":
		with open(workpath+"/Shellscript/subgroupdir2.sh", "w") as subdirfile:
			subdirfile.write("#!/usr/bin/sh\n"+groupsh)
			subdirfile.close()
			if not os.path.exists(workpath+"/07_groups/group_"+"_vs_".join(groupes)) and os.path.exists(workpath+"/07_groups/"):
				os.system("sh "+workpath+"/Shellscript/subgroupdir2.sh")

def samplefastq(samplefile):
	work_path = os.getcwd()
	sample2fastq = {}
	with open(samplefile, 'rt') as fin:
		linenum = 0
		for line in fin:
			linenum += 1
			if linenum > 1:
				tmp = line.replace("\n","").split("\t")
				sample = tmp[0]
				files = tmp[1].split(";")
				sample2fastq[sample] = files
	return sample2fastq

def congfile(configurefile):
	filedict = {}
	with open(configurefile, 'rt') as fin:
		for line in fin:
			tmp = line.replace("\n", "").split("=")
			filedict[tmp[0]] = tmp[1]
	return filedict

def quality_check(workpath, quality_out, readslist, fqrun):
	raw_clean = quality_out
	fqsh = []
	for i in xrange(len(readslist)):
		fqsh.append("fastqc -q --extract -o fastqc_out "+readslist[i])
	fqtmp = "\n".join(fqsh)
	if raw_clean == "raw":
		if not os.path.exists(workpath+"/00_raw_reads/fastqc_out_run.log"):	
			with open(workpath+"/Shellscript/00_raw_reads.sh" , "w") as shfile:
				shfile.write("#!/usr/bin/sh\ncd "+workpath+"/00_raw_reads\n" + fqtmp+"\ntouch fastqc_out_run.log")
				shfile.close()
				if fqrun == "Yes":
					os.system("sh "+workpath+"/Shellscript/00_raw_reads.sh")
	if raw_clean == "clean":
		if not os.path.exists(workpath+"/01_clean_reads/fastqc_out_run.log"):
			with open(workpath+"/Shellscript/02_clean_reads.sh" , "w") as shfile:
				shfile.write("#!/usr/bin/sh\ncd "+workpath+"/01_clean_reads\n" + fqtmp+"\ntouch fastqc_out_run.log")
				shfile.close()
				if fqrun == "Yes":
					os.system("sh "+workpath+"/Shellscript/02_clean_reads.sh")

def pick_quality(workpath, samplelist):
	copysh = []
	for sample in samplelist:
		copyimage = "cp "+workpath+"/01_clean_reads/fastqc_out/"+sample+"_clean_merge_fastqc/Images/per_base_quality.png "+workpath+"/01_clean_reads/quality_image/"+sample+"_quality.png"
		copysh.append(copyimage)
	if not os.path.exists(workpath+"01_clean_reads/01_pick_quality.done") and os.path.exists(workpath+"01_clean_reads/fastqc_out"):
		with open(workpath+"/Shellscript/02_clean_quality_pick.sh" , "w") as shfile:
			shfile.write("#!/usr/bin/sh\ncd "+workpath+"/01_clean_reads\nmkdir quality_image \n"+"\n".join(copysh))
			shfile.close()
			os.system("sh "+workpath+"/Shellscript/02_clean_quality_pick.sh")

			
def assembly_raw(workpath, Adapter_Primer3, Adapter_Primer5 , sample2readsfile, maxoverlength, mismatch, cleandatapath, minlength, run_Assemble):
	if not os.path.exists(workpath+"/02_assemble/01_raw_merge.log") and os.path.exists(workpath+"/00_raw_reads/fastqc_out"):
		assembly = []
		cutadaptcmd = []
		for sample in sample2readsfile:
			cutadaptstr = "cutadapt -a "+ Adapter_Primer3 +" "+cleandatapath+"/"+sample+"_raw_merge.fastq --quiet -o "+cleandatapath + "/" + sample +"_clean_merge_tmp.fastq\ncutadapt -g "+ Adapter_Primer5+" "+cleandatapath+"/"+sample +"_clean_merge_tmp.fastq --max-n 2 -m " + str(minlength) + " -o "+cleandatapath+"/"+sample+"_clean_merge.fastq\nwc "+cleandatapath+"/"+sample+"_clean_merge.fastq | awk -v tmp='"+sample+"' '{print tmp\"\\t\"$1/4}' >> Clean_Nums.txt" 
			cutadaptcmd.append(cutadaptstr)
			if len(sample2readsfile[sample]) == 2:
				tmpstr = "flash -q -c -M "+str(maxoverlength)+" -x "+str(mismatch)
				readstmp = sample2readsfile[sample][1] + " " + sample2readsfile[sample][0] +" > "+cleandatapath+"/"+sample+"_raw_merge.fastq"
				Rawnum = "wc "+sample2readsfile[sample][0]+" | awk -v tmp='"+sample+"' '{print tmp\"\t\"$1/4}' >> Raw_Nums.txt"
				assembly.append(tmpstr+" "+readstmp+"\n"+Rawnum)
			if len(sample2readsfile[sample]) == 1:
				tmpstr = "cp "+sample2readsfile[sample][0]+" "+cleandatapath+"/"+sample+"_raw_merge.fastq\nwc "+cleandatapath+"/"+sample+"_raw_merge.fastq | awk -v tmp='"+sample+"' '{print tmp\"\\t\"$1/4}' >> Raw_Nums.txt"
				assembly.append(tmpstr)
		with open(workpath+"/Shellscript/01_assemble.sh" , "w") as shfile:
			shfile.write("#!/usr/bin/sh\ncd "+workpath+"/02_assemble/\n"+"\n".join(assembly)+"\n"+"\n".join(cutadaptcmd)+"\ntouch 01_raw_merge.log")	
			shfile.close()
			if run_Assemble == "Yes":
				os.system("sh "+workpath+"/Shellscript/01_assemble.sh")

def totalsamplemerge(workpath, sample2readsfile, cleandatapath, run_merge_rename):
	if not os.path.exists(workpath+"/02_assemble/02_total_merge.log") and os.path.exists(workpath+"/01_clean_reads/fastqc_out"):
		mergecat = []
		qiimecat = []
		for sample in sample2readsfile:
			awkstr = "awk -v tmp=\""+sample+"\" 'BEGIN{n=1;j=0}{if(n%4==1){j++;print \">\"tmp\".\"j\";barcodelabel=\"tmp}else{if(n%4==2){print $0}};n++}' "+cleandatapath+"/"+sample+"_clean_merge.fastq" + " >> 02_assemble/usearch_input_merge.fa"
			qiimestr = "awk -v tmp=\""+sample+"\" 'BEGIN{n=1;j=0}{if(n%4==1){j++;print \">\"tmp\"_\"j}else{if(n%4==2){print $0}};n++}' "+cleandatapath+"/"+sample+"_clean_merge.fastq" + " >> 02_assemble/qiime_input_merge.fa"
			mergecat.append(awkstr)
			qiimecat.append(qiimestr)
		with open(workpath+"/Shellscript/02_otu_input.sh" , "w") as shfile:
			shfile.write("#!/usr/bin/sh\n"+"\n".join(mergecat)+"\n"+"\n".join(qiimecat)+"\ntouch 02_assemble/02_total_merge.log")
			shfile.close()
			if run_merge_rename == "Yes":
				os.system("sh "+workpath+"/Shellscript/02_otu_input.sh")

def pick_otu_my(workpath, similarity, gold_database, runPickOTU):
	if not os.path.exists(workpath+"/03_pick_OTU/03_pick_OTU.log") and os.path.exists(workpath+"/02_assemble/usearch_input_merge.fa"):
		with open(workpath+"/Shellscript/03_pick_otu.sh" , "w") as shfile:
			shfile.write("#!/usr/bin/sh\ncd "+workpath+"/03_pick_OTU/\nusearch10 -fastx_uniques "+workpath+"/02_assemble/usearch_input_merge.fa --fastaout derep_mix_clean.fa --sizeout\nusearch10 -sortbysize derep_mix_clean.fa -fastaout derep_mix_clean.sorted.fa -minsize 2\nusearch10 -cluster_otus derep_mix_clean.sorted.fa -otus otus1.fa -relabel OTU_ -uparseout results.txt\nusearch8 -uchime_ref otus1.fa -db "+gold_database+" -strand plus -nonchimeras otus.fa\nusearch10 -usearch_global "+workpath+"/02_assemble/usearch_input_merge.fa -db otus1.fa -strand plus -id "+ str(similarity) +" -uc reads_map1.uc\nusearch10 -usearch_global "+workpath+"/02_assemble/usearch_input_merge.fa -db otus.fa -strand plus -id " + str(similarity) +" -uc reads_map.uc\npython ~/software/16s/usearch/python_scripts/uc2otutab.py reads_map1.uc > otus_table1.txt\npython ~/software/16s/usearch/python_scripts/uc2otutab.py reads_map.uc > otus_table.txt\ntouch 03_pick_OTU.log\nBacteria_assemble_stastic.R "+workpath+"/02_assemble/Raw_Nums.txt "+workpath+"/02_assemble/Clean_Nums.txt otus_table1.txt otus_table.txt")
			shfile.close()
			if runPickOTU == "Yes":
				os.system("sh "+workpath+"/Shellscript/03_pick_otu.sh")

def taxannotation(workpath, taxonomy, tax_fasta, projecttype, rdpmemory):
	#if not os.path.exists(workpath+"/04_OTU_assignment/04_tax_annotaion.log") and os.path.exists(workpath+"/03_pick_OTU/otus.fa"):
	with open(workpath+"/Shellscript/04_otu_tax_annotation.sh", "w") as shfile:
		if projecttype=="16s" or projecttype=="16S" or projecttype=="18S" or projecttype=="18s":
			shfile.write("#!/usr/bin/sh\ncd "+workpath+"/04_OTU_assignment\ncp "+workpath+"/03_pick_OTU/otus.fa .\nassign_taxonomy.py -i otus.fa -t "+ taxonomy + " -r " + tax_fasta + " -o otus_assign \ntouch 04_tax_annotaion.log")
			shfile.close()
		elif projecttype=="its" or projecttype=="ITS":
			shfile.write("#!/usr/bin/sh\ncd "+workpath+"/04_OTU_assignment\ncp "+workpath+"/03_pick_OTU/otus.fa .\nassign_taxonomy.py -m rdp -i otus.fa -t "+ taxonomy + " -r " + tax_fasta + " --rdp_max_memory "+str(rdpmemory)+" -o otus_assign \ntouch 04_tax_annotaion.log")	
	if not os.path.exists(workpath+"/04_OTU_assignment/04_tax_annotaion.log") and os.path.exists(workpath+"/03_pick_OTU/otus.fa"):
		os.system("sh "+workpath+"/Shellscript/04_otu_tax_annotation.sh")

def OTUtabBuild(workpath, Kindom):
	taxlevel = ["Phylum", "Class", "Order", "Family", "Genus"]
	if not os.path.exists(workpath+"/05_OTU_table_build/05_OTU_tab.log") and os.path.exists(workpath+"/04_OTU_assignment/otus_assign"):
		cmdtmplist = []
		cmdtmpstr = ""
		for i in xrange(len(taxlevel)):
			l = i+2
			cmd = "taxonomy_class_count.R otu_tax_table.txt -o "+taxlevel[i]+"_data.xls -l "+str(l)
			cmdtmplist.append(cmd)
		cmdtmpstr = "\n".join(cmdtmplist)
		with open(workpath+"/Shellscript/05_OTUtab_Build.sh", "w") as shfile:			
			if Kindom == "Bacteria" or Kindom == "bacteria" or Kindom=="":
				shfile.write("#!/usr/bin/sh\ncd " + workpath + "/05_OTU_table_build\notu2reads.py "+workpath+"/03_pick_OTU/reads_map.uc\ncp "+workpath+"/04_OTU_assignment/otus_assign/otus_tax_assignments.txt .\nmake_otu_table.py -i otu2seqs.tab -t otus_tax_assignments.txt -o otu_tax_table.biom\nbiom convert -i otu_tax_table.biom -o otu_tax_table.txt --to-tsv --header-key=taxonomy\nqiime_close_otutax.py otu_tax_table.txt\ngrep -v 'k__Unassigned;' otu_tax_table.tab > otu_tax_table.txt\nrare_otu_count.R otu_tax_table.txt -o otu_tax_table.tab\nmv otu_tax_table.tab otu_tax_table.txt\nrm otu_tax_table.biom\ntxt2biom.py otu_tax_table.txt otu_tax_table.biom\nrare_norm.R otu_tax_table.txt -o otu_tax_table_norm.txt \ntxt2biom.py otu_tax_table_norm.txt otu_tax_table_norm.biom\nextract_seq.py "+workpath+"/03_pick_OTU/otus.fa otu_tax_table.txt > ../03_pick_OTU/rep_seq.fa\ntouch 05_OTU_tab.log\n"+cmdtmpstr+"\ncd "+workpath+"/03_pick_OTU/\nalign_seqs.py -m muscle -i rep_seq.fa -o ./\nmake_phylogeny.py -i rep_seq_aligned.fasta -o rep_seq_phylo.tre")
				shfile.close()
			elif Kindom == "Archaea" or Kindom == "archaea":
				shfile.write("#!/usr/bin/sh\ncd " + workpath + "/05_OTU_table_build\notu2reads.py "+workpath+"/03_pick_OTU/reads_map.uc\ncp "+workpath+"/04_OTU_assignment/otus_assign/otus_tax_assignments.txt .\nmake_otu_table.py -i otu2seqs.tab -t otus_tax_assignments.txt -o otu_tax_table.biom\nbiom convert -i otu_tax_table.biom -o otu_tax_table.txt --to-tsv --header-key=taxonomy\nqiime_close_otutax.py otu_tax_table.txt\ngrep -v 'k__Unassigned;' otu_tax_table.tab | grep -v 'k__Bacteria;' > otu_tax_table.txt\nrm otu_tax_table.biom\ntxt2biom.py otu_tax_table.txt otu_tax_table.biom\n"+cmdtmpstr+"\nArchaea_assemble_stastic.R "+workpath+"/02_assemble/Raw_Nums.txt "+workpath+"/02_assemble/Clean_Nums.txt "+workpath+"/03_pick_OTU/otus_table1.txt "+workpath+"/03_pick_OTU/otus_table.txt otu_tax_table.txt\nmv Statistic_RawCleanNums.xls "+workpath+"/03_pick_OTU/\nrare_otu_count.R otu_tax_table.txt -o otu_tax_table.tab\nmv otu_tax_table.tab otu_tax_table.txt\nrare_norm.R otu_tax_table.txt -o otu_tax_table_norm.txt\ntxt2biom.py otu_tax_table_norm.txt otu_tax_table_norm.biom\nextract_seq.py "+workpath+"/03_pick_OTU/otus.fa otu_tax_table.txt > ../03_pick_OTU/rep_seq.fa\ntouch 05_OTU_tab.log\n"+cmdtmpstr+"\ncd "+workpath+"/03_pick_OTU/\nalign_seqs.py -m muscle -i rep_seq.fa -o ./\nmake_phylogeny.py -i rep_seq_aligned.fasta -o rep_seq_phylo.tre")
				shfile.close()
			os.system("sh "+workpath+"/Shellscript/05_OTUtab_Build.sh")

def rare_analysis(workpath, sample_qiime):
	if not os.path.exists(workpath+"/06_rare_analysis/06_rare_analysis.log") and os.path.exists(workpath+"/05_OTU_table_build/otu_tax_table.biom"):
		with open(workpath+"/Shellscript/06_rare_analysis.sh", "w") as shfile:
			shfile.write("#!/usr/bin/sh\ncd "+workpath +"/06_rare_analysis/\necho -e \"alpha_diversity:metrics\tshannon,chao1,simpson,goods_coverage,observed_species\" > alpha_params.txt\ncp "+workpath+"/05_OTU_table_build/otu_tax_table.biom .\nalpha_rarefaction.py -i otu_tax_table.biom -f -p alpha_params.txt -n 10 -o arare -m "+sample_qiime+"\nfor file in `ls arare/alpha_div_collated/*txt`; do tmp1=${file/.txt/};tmp=${tmp1/arare\/alpha_div_collated\//};rarefraction_smooth.R $file $tmp $tmp;done\ntouch 06_rare_analysis.log")
			shfile.close()
			os.system("sh "+workpath+"/Shellscript/06_rare_analysis.sh")

def subgroup(workpath, groups, sample_total, runsplit, group):
	tmpvs = "_vs_".join(groups.split(","))
	subgroup = "subsample_"+tmpvs+".txt"
	flagfile = workpath+"/Tmpgroups/group_"+tmpvs+"/"+subgroup
	flagrun = runsplit
	if os.path.exists(workpath+"/05_OTU_table_build/otu_tax_table.biom") and not os.path.exists(flagfile):
		tmpsh = workpath+"/Shellscript/07_subgroup"+tmpvs+".sh"
		if group== "group":
			with open(tmpsh, "w") as shfile:
				shfile.write("#!/usr/bin/sh\nsubsample_file.R "+sample_total+" "+groups+" -o "+flagfile)
				shfile.close()
				if flagrun== "Yes": 
					os.system("sh "+tmpsh)
		if group=="subgroup":
			with open(tmpsh, "w") as shfile:
				shfile.write("#!/usr/bin/sh\nsubsample_split.R " + sample_total + " " +groups+" -o "+flagfile)
				shfile.close()
				if flagrun== "Yes":
					os.system("sh "+tmpsh)

def splitOTUtab(workpath, groups, subanalysis):
	tmpvs = "_vs_".join(groups.split(","))
	subTab = "otu_tax_tab_"+tmpvs+".biom"
	subTabnorm = "otu_tab_norm_"+tmpvs+".biom"
	subtxt = "otu_tax_tab_"+tmpvs+".txt"
	subtxtnorm = "otu_tab_norm_"+tmpvs+".txt"
	flagfile =  workpath+"/07_groups/group_"+tmpvs+"/"+subTab
	flagfilenorm = workpath+"/07_groups/group_"+tmpvs+"/"+subTabnorm
	filetxt = workpath+"/07_groups/group_"+tmpvs+"/"+subtxt
	filetxtnorm = workpath+"/07_groups/group_"+tmpvs+"/"+subtxtnorm
	tmpsamp = workpath+"/Tmpgroups/group_"+tmpvs+"/subsample_"+tmpvs+".txt"
	if os.path.exists("05_OTU_table_build/otu_tax_table.biom") and not os.path.exists(flagfile):
		tmpsh = workpath+"/Tmpgroups/group_"+tmpvs+"/split_tab"+tmpvs+".sh"
		with open(tmpsh, "w") as shfile:
			shfile.write("#!/usr/bin/sh\nfilter_samples_from_otu_table.py -i "+workpath+"/05_OTU_table_build/otu_tax_table.biom -o "+flagfile+" --sample_id_fp "+tmpsamp+"\nbiom convert -i "+flagfile+" -o "+filetxt+" --header-key=taxonomy --to-tsv\nrare_norm.R "+filetxt+" -o "+filetxtnorm+"\ntxt2biom.py "+filetxtnorm+" "+flagfilenorm)
			shfile.close()
			if subanalysis == "Yes":
				os.system("sh "+tmpsh)

def Vennplot(workpath, groups, subanalysis):
	tmpvs = "_vs_".join(groups.split(","))
	Vennimage = "Venn_plot_"+tmpvs+".svg"
	Anoimage = "Anosim_plot_"+tmpvs+".svg"
	subtxtnorm = "otu_tab_norm_"+tmpvs+".txt"
	tmpsam = workpath+"/Tmpgroups/group_"+tmpvs+"/subsample_"+tmpvs+".txt"
	filetxtnorm = workpath+"/07_groups/group_"+tmpvs+"/"+subtxtnorm
	outimage = workpath+"/07_groups/group_"+tmpvs+"/"+Vennimage
	Anosim = workpath+"/07_groups/group_"+tmpvs+"/08_Anosim_analysis/"+ Anoimage
	if os.path.exists(filetxtnorm) and not os.path.exists(Anosim):
		tmpsh = workpath+"/Tmpgroups/group_"+tmpvs+"/08_Anosim_Venn"+tmpvs+".sh"
		with open(tmpsh, "w") as shfile:
			shfile.write("#!/usr/bin/sh\nvenn_plot_plus.R "+filetxtnorm+" "+tmpsam+" -o "+outimage+"\nAnosim_plot.R "+filetxtnorm+" "+tmpsam+" -o "+Anosim)
			shfile.close()
			if subanalysis =="Yes":
				os.system("sh "+tmpsh)

def alpha_plot(workpath, groups, subanalysis):
	tmpvs = "_vs_".join(groups.split(","))
	alphaimage = "alpha_plot_"+tmpvs+".svg"
	alphatab = "alpha_tab_" + tmpvs +".txt"
	subtxt = "otu_tax_tab_"+tmpvs+".txt"
	tmpsam = workpath+"/Tmpgroups/group_"+tmpvs+"/subsample_"+tmpvs+".txt"
	filetxt = workpath+"/07_groups/group_"+tmpvs+"/"+subtxt
	outimage = workpath+"/07_groups/group_"+tmpvs+"/09_alpha_analysis/"+alphaimage
	outdata = workpath+"/07_groups/group_"+tmpvs+"/09_alpha_analysis/"+ alphatab
	if os.path.exists(filetxt) and not os.path.exists(outimage):
		tmpsh = workpath+"/Tmpgroups/group_"+tmpvs+"/09_alpha_plot"+tmpvs+".sh"
		with open(tmpsh, "w") as shfile:
			shfile.write("#!/usr/bin/sh\nalpha_rare_plot.R "+filetxt+" "+tmpsam+" -d "+ outdata+" -o "+outimage)
			shfile.close()
			if subanalysis == "Yes":
				os.system("sh "+tmpsh)

def pca_plot(workpath, groups, subanalysis):
	tmpvs = "_vs_".join(groups.split(","))
	pcaimage = "pca_plot_"+tmpvs+".svg"
	subtxtnorm = "otu_tab_norm_"+tmpvs+".txt"
	tmpsam = workpath+"/Tmpgroups/group_"+tmpvs+"/subsample_"+tmpvs+".txt"
	outimage = workpath+"/07_groups/group_"+tmpvs+"/10_beta_analysis/"+pcaimage
	filetxtnorm = workpath+"/07_groups/group_"+tmpvs+"/"+subtxtnorm
	if os.path.exists(filetxtnorm) and not os.path.exists(outimage):
		tmpsh = workpath+"/Tmpgroups/group_"+tmpvs+"/10_pca_plot"+tmpvs+".sh"
		with open(tmpsh, "w") as shfile:
			shfile.write("#!/usr/bin/sh\npca_box_plot.R "+filetxtnorm+" "+tmpsam+" -o "+outimage)
			shfile.close()
			if subanalysis == "Yes":
				os.system("sh "+tmpsh)

def PCoA_NMDS_plot(workpath, groups, subanalysis):
	tmpvs = "_vs_".join(groups.split(","))
	PCoAimage = "PCoA_plot_" + tmpvs + ".svg"
	NMDSimage = "NMDS_plot_" + tmpvs + ".svg"
	subbiomnorm = "otu_tab_norm_" + tmpvs + ".biom"
	tmpsam = workpath + "/Tmpgroups/group_" + tmpvs + "/subsample_" + tmpvs + ".txt"
	PCoAoutimage = workpath + "/07_groups/group_" + tmpvs + "/10_beta_analysis/" + PCoAimage
	NMDSoutimage = workpath + "/07_groups/group_" + tmpvs + "/10_beta_analysis/" + NMDSimage
	fileBiomnorm = workpath+"/07_groups/group_"+tmpvs+"/"+subbiomnorm
	if os.path.exists(fileBiomnorm) and not os.path.exists(PCoAoutimage):
		tmpsh = workpath+"/Tmpgroups/group_"+tmpvs+"/10_PCoA_NMDS_plot"+tmpvs+".sh"
		with open(tmpsh, "w") as shfile:
			shfile.write("#!/usr/bin/sh\nPCoA_NMDS_phyloseq_plot.R " + fileBiomnorm + " " + tmpsam + " -s 2 -m2 " + PCoAoutimage + " -m3 " + NMDSoutimage)
			shfile.close()
			if subanalysis == "Yes":
				os.system("sh "+tmpsh)

def taxclass(workpath, groups, subanalysis, heatmap_width, heatmap_height):
	tmpvs = "_vs_".join(groups.split(","))
	subtxtnorm = "otu_tab_norm_"+tmpvs+".txt"
	filetxtnorm = workpath+"/07_groups/group_"+tmpvs+"/"+subtxtnorm
	tmpsam = workpath+"/Tmpgroups/group_"+tmpvs+"/subsample_"+tmpvs+".txt"
	tmpshlist = []
	htshlist = []
	htshlistno = []
	for i in xrange(2,8):
		taximage = "otu_taxlevel_"+str(i)+".svg"
		heatmapeout = "heatmap_level_"+str(i)+".svg"
		heatmapeoutno = "heatmap_level_"+str(i)+"_noclust.svg"
		taxtxt = "otu_taxlevel_"+str(i)+".txt"
		outimage = workpath+"/07_groups/group_"+tmpvs+"/11_tax_classification/"+taximage
		outdata = workpath+"/07_groups/group_"+tmpvs+"/11_tax_classification/"+taxtxt
		heatmapimage = workpath+"/07_groups/group_"+tmpvs+"/11_tax_classification/"+heatmapeout
		heatmapimageno = workpath+"/07_groups/group_"+tmpvs+"/11_tax_classification/"+heatmapeoutno
		tmpshlist.append("taxonomy_bar_plot.R "+filetxtnorm+" "+tmpsam+" -l "+str(i)+" -d "+outdata+" -o "+outimage)
		htshlist.append("heatmap_plot.R "+outdata+" "+tmpsam+" -w "+str(heatmap_width)+" -H "+str(heatmap_height)+" -T N -o "+ heatmapimage)
		htshlistno.append("heatmap_plot.R "+outdata+" "+tmpsam+" -w "+str(heatmap_width)+" -H "+str(heatmap_height)+" -T N -F N -o "+ heatmapimageno)
	if os.path.exists(filetxtnorm) and not os.path.exists(outimage):
		tmpsh = workpath+"/Tmpgroups/group_"+tmpvs+"/11_tax_class"+tmpvs+".sh"
		with open(tmpsh, "w") as shfile:
			shfile.write("#!/usr/bin/sh\n"+"\n".join(tmpshlist)+"\n"+"\n".join(htshlist)+"\n"+"\n".join(htshlistno))
			shfile.close()
			if subanalysis == "Yes":
				os.system("sh "+tmpsh)

def cluster_analysis(workpath, groups, subanalysis, cl_label_size, clbar_tree_height, blank_width, clbar_labelsize):
	tmpvs = "_vs_".join(groups.split(","))
	subtxtnorm = "otu_tab_norm_"+tmpvs+".txt"
	filetxtnorm = workpath+"/07_groups/group_"+tmpvs+"/"+subtxtnorm
	tmpsam = workpath+"/Tmpgroups/group_"+tmpvs+"/subsample_"+tmpvs+".txt"
	clustdone = workpath+"/Tmpgroups/group_"+tmpvs+"/cluster_run.done"
	clshlist = []
	clbarshlist = []
	for i in xrange(2,8):
		climage = "otu_taxlevel_"+str(i)+"_clust.svg"
		clbarimage = "otu_taxlevel_"+str(i)+"_clust_bar.svg"
		outimage = workpath+"/07_groups/group_"+tmpvs+"/11_tax_classification/"+climage
		clbaroutimage = workpath+"/07_groups/group_"+tmpvs+"/11_tax_classification/"+clbarimage
		clshlist.append("cluster_plot.R "+filetxtnorm+" "+tmpsam+" -l "+str(i)+" -s "+str(cl_label_size)+" -o "+outimage)
		clbarshlist.append("cluster_plot_plus.R "+filetxtnorm+" "+tmpsam+" -l "+str(i)+" -a "+str(clbar_tree_height)+" -z "+str(blank_width)+" -s "+str(clbar_labelsize)+" -o "+clbaroutimage)
	if os.path.exists(filetxtnorm) and not os.path.exists(clustdone):
		tmpsh = workpath+"/Tmpgroups/group_"+tmpvs+"/11_clust_run_"+tmpvs+".sh"
		with open(tmpsh, "w") as shfile:
			shfile.write("#!/usr/bin/sh\n"+"\n".join(clshlist)+"\n"+"\n".join(clbarshlist)+"\ntouch "+clustdone)
			shfile.close()
			if subanalysis == "Yes":
				os.system("sh "+tmpsh)

def LEfSe_DOTU(workpath, groups, subanalysis):
	tmpvs = "_vs_".join(groups.split(","))
	subTabnorm = "otu_tab_norm_"+tmpvs+".biom"
	flagfilenorm = workpath+"/07_groups/group_"+tmpvs+"/"+subTabnorm
	tmpsam = workpath+"/Tmpgroups/group_"+tmpvs+"/subsample_"+tmpvs+".txt"
	outpath = workpath+"/07_groups/group_"+tmpvs+"/12_LEfSe_analysis/"
	tmppath = workpath+"/Tmpgroups/group_"+tmpvs+"/"
	LEfSe_input = outpath+"LEfSe_input_"+tmpvs+".txt"
	LEfSe_tmp_input = outpath+"LEfSe_tmp_input_"+tmpvs+".txt"
	LEfSe_input_format = tmppath+"LEfSe_input_"+tmpvs+"_format.txt"
	LEfSe_formatres = tmppath+"LEfSe_input_"+tmpvs+"_format.txt.res"
	LEfSe_DEtaxtab = outpath+"DEtax_LEfSe_tab"+tmpvs+".xls"
	LEfSe_LDA = outpath+"LEfSe_LDA_plot_"+tmpvs+".svg"
	LEfSe_caladogram = outpath+"LEfSe_caladogram_"+tmpvs+".svg"
	biomaker = outpath+"/biomarkers_"+tmpvs
	if os.path.exists(flagfilenorm) and not os.path.exists(LEfSe_caladogram):
		tmpsh = workpath+"/Tmpgroups/group_"+tmpvs+"/12_LEfSe_analysis_"+tmpvs+".sh"
		with open(tmpsh, 'w') as shfile:
			shfile.write("#!/usr/bin/sh\nQiimeToLEfSe.R "+flagfilenorm+" "+tmpsam+" -d "+LEfSe_tmp_input+" -o "+LEfSe_input+"\nformat_input.py "+LEfSe_input+" "+LEfSe_input_format+" -c 1 -o 1000000\nrun_lefse.py "+LEfSe_input_format+" "+LEfSe_formatres+"\nLDA_cicular_plot.R "+LEfSe_formatres+" "+tmpsam+" -o "+LEfSe_LDA+"\nplot_cladogram.py "+LEfSe_formatres+" "+LEfSe_caladogram+" --format svg --dpi 400 --left_space_prop 0.15 --right_space_prop 0.3\nmkdir -p "+biomaker+"\nplot_features.py "+LEfSe_input_format+" "+LEfSe_formatres+" "+biomaker+"/ --title_font_size 10 --dpi 300\nLEfSe_DEtax_tab.R "+LEfSe_formatres+" "+LEfSe_tmp_input+" -o "+LEfSe_DEtaxtab+"\ncd "+biomaker+"\ni=1;for file in `ls *.png`;do mv $file 'LEfSe_DET_Figure'$i'.png';i=$[$i+1];done")
			shfile.close()
			if subanalysis == "Yes":
				os.system("sh "+tmpsh)

def EnvironRDA(workpath, groups):
	tmpvs = "_vs_".join(groups.split(","))
	subtxtnorm = "otu_tab_norm_"+tmpvs+".txt"
	filetxtnorm = workpath+"/07_groups/group_"+tmpvs+"/"+subtxtnorm
	tmpsam = workpath+"/Tmpgroups/group_"+tmpvs+"/subsample_"+tmpvs+".txt"
	RDAimage = "RDA_plot_"+tmpvs+".svg"
	CCAimage = "CCA_plot_"+tmpvs+".svg"
	outpath = workpath+"/07_groups/group_"+tmpvs+"/12_LEfSe_analysis/"
	RDAoutimage = outpath+RDAimage
	CCAoutimage = outpath+CCAimage
	if os.path.exists(filetxtnorm) and not os.path.exists(RDAoutimage):
		tmpsh = workpath+"/Tmpgroups/group_"+tmpvs+"/13_enrironment_analysis_"+tmpvs+".sh"
		with open(tmpsh, 'w') as shfile:
			shfile.write("#!/usr/bin/sh\nRscript RDA_CCA_environment.R "+filetxtnorm+" "+tmpsam+" environment_factor.txt -r "+RDAoutimage+" -c "+CCAoutimage)
			shfile.close()
			
def myspeciestree(workpath, groups, subanalysis):
	tmpvs = "_vs_".join(groups.split(","))
	subTabnorm = "otu_tab_norm_"+tmpvs+".biom"
	flagfilenorm = workpath+"/07_groups/group_"+tmpvs+"/"+subTabnorm
	tmpsam = workpath+"/Tmpgroups/group_"+tmpvs+"/subsample_"+tmpvs+".txt"
	outpath = workpath+"/07_groups/group_"+tmpvs+"/14_species_tree_distribution/"
	tmppath = workpath+"/Tmpgroups/group_"+tmpvs+"/"
	inputfile = outpath+"Tree_input_"+tmpvs+".txt"
	tmpfile = tmppath+"Tree_input_tmp1_"+tmpvs+".txt"
	treefile = outpath+"Tree_input_"+tmpvs+".tree"
	tmpannot = tmppath+"tmp.annot"
	annotfile = outpath+"Tree_input_"+tmpvs+".annot"
	xmlannot = outpath+"Tree_input_"+tmpvs+".xml"
	treeoutimage = outpath+"Species_Tree_"+tmpvs+".svg"
	if os.path.exists(flagfilenorm) and not os.path.exists(treeoutimage):
		tmpsh = workpath+"/Tmpgroups/group_"+tmpvs+"/14_Tree_plot_"+tmpvs+".sh"
		with open(tmpsh, 'w') as shfile:
			shfile.write("#!/usr/bin/sh\nQiimeToLEfSe.R "+flagfilenorm+" "+tmpsam+" -o "+inputfile+"\npick_known_tax_to_tree.py "+inputfile+" g > "+tmpfile+"\nmetaphlan2graphlan.py "+tmpfile+" --tree_file "+treefile+" --annot_file "+tmpannot+"\ngenerate_annotaion_for_species_tree.py "+tmpfile+" > "+annotfile+"\ngraphlan_annotate.py "+treefile+" "+xmlannot+" --annot "+annotfile+"\ngraphlan.py "+xmlannot+" "+treeoutimage+" --size 16")
			shfile.close()
			if subanalysis == "Yes":
				os.system("sh "+tmpsh)

def PDFbuildtotal(workpath):
	#tmpvs = "_vs_".join(groups.split(","))
	pickegglog = workpath+"/08_Function_predict/08_pick_close_OTU.log"
	pdfdonelog = workpath+"/09PDFbuild/09PDFdone1.log"
	if os.path.exists(pickegglog) and not os.path.exists(pdfdonelog):
		tmpsh = workpath+"/Shellscript/09PDFbuild01.sh"	
		with open(tmpsh, "w") as shfile:
			shfile.write("#!/usr/bin/sh\nPDF_build_16S.sh\ntouch 09PDFbuild/09PDFdone1.log")
			shfile.close()
			os.system("sh "+tmpsh)

def PDFbuildgroup(workpath, groups):
	tmpvs = "_vs_".join(groups.split(","))
	pickegglog = workpath+"/08_Function_predict/08_pick_close_OTU.log"
	pdfdone2log = workpath+"/09PDFbuild/09PDFdone2.log"
	if os.path.exists(pickegglog):#and not os.path.exists(pdfdone2log):
		tmpsh = workpath+"/Shellscript/09PDFbuild02.sh"
		with open(tmpsh, "w") as shfile:
			shfile.write("#!/usr/bin/sh\nselect_group_floder_16S.sh "+tmpvs)
			shfile.close()
			os.system("sh "+tmpsh)		


if __name__ == "__main__":
	params = readparas(sys.argv)
	sampleinfo = params['sample_file']
	configure = params['configure_file']
	confiles = congfile(configure)
	workpath = confiles['workpath']
	datapath = confiles['datapath']
	samplefile = samplefastq(sampleinfo)
	projecttype = confiles["projecttype"]
	Kindom = confiles["Kindom"]
	Adapter_Primer3 = confiles["Adapter_Primer3"]
	Adapter_Primer5 = confiles['Adapter_Primer5']
	rdpmemory = confiles["rdpmemory"]
	QCrun = confiles['QCrun']
	groups = confiles["groups"]
	subgroups = confiles["subgroups"]
	maxoverlength = confiles["maxoverlength"]
	mismatch = confiles["mismatch"]
	cleandatapath = confiles["cleandatapath"]
	minlength = confiles['minlength']
	run_Assemble = confiles['run_Assemble']
	run_merge_rename = confiles['run_merge_rename']
	gold_database = confiles['gold_database']
	similarity = confiles["similarity"]
	runPickOTU = confiles["runPickOTU"]
	taxonomy = confiles['taxonomy']
	tax_fasta = confiles['tax_fasta']
	sample_qiime = confiles['sample_qiime']
	heatmap_width = confiles['heatmap_width']
	heatmap_height = confiles['heatmap_height']
	cl_label_size = confiles["cl_label_size"]
	clbar_tree_height = confiles["clbar_tree_height"]
	blank_width = confiles["blank_width"]
	clbar_labelsize = confiles["clbar_labelsize"]
	environment = confiles["environmentfactor"]
	runsplit = confiles["runsplit"]
	subanalysis = confiles["subanalysis"]
	mkdir_16s()
	#print len(groups.split(";")[0])
	if len(groups.split(";")[0]) !=0 :
		subgroupdir(workpath, groups, "group")
	if len(subgroups.split(";")[0]) != 0:
		subgroupdir(workpath, subgroups, "subgroup")
	readslist = []
	sample2readsfile = {}
	samplelist = samplefile.keys()
	for sample in samplefile:
		if len(samplefile[sample]) == 2:
			readsfile = datapath+"/"+samplefile[sample][0]+" "+datapath+"/"+samplefile[sample][1]
			sample2readsfile[sample] = []
			sample2readsfile[sample].append(datapath+"/"+samplefile[sample][0])
			sample2readsfile[sample].append(datapath+"/"+samplefile[sample][1])
		if len(samplefile[sample]) == 1:
			readsfile = datapath+"/"+samplefile[sample][0]
			sample2readsfile[sample] = []
			sample2readsfile[sample].append(datapath+"/"+samplefile[sample][0])
		readslist.append(readsfile)
	cleanreadslist = []
	for sample in sample2readsfile:
		cleanreadsfile = cleandatapath+"/"+sample+"_clean_merge.fastq"
		cleanreadslist.append(cleanreadsfile)
	quality_check(workpath, "raw", readslist, QCrun)
	assembly_raw(workpath, Adapter_Primer3, Adapter_Primer5, sample2readsfile, maxoverlength, mismatch, cleandatapath, minlength, run_Assemble)
	quality_check(workpath, "clean", cleanreadslist, QCrun)
	pick_quality(workpath, samplelist)
	totalsamplemerge(workpath, sample2readsfile, cleandatapath, run_merge_rename)
	pick_otu_my(workpath, similarity, gold_database, runPickOTU)
	taxannotation(workpath, taxonomy, tax_fasta, projecttype, rdpmemory)
	OTUtabBuild(workpath, Kindom)
	rare_analysis(workpath, sample_qiime)
	PDFbuildtotal(workpath)
	if len(groups.split(";")[0]) != 0:
		groupslist = groups.split(";")
		for i in xrange(len(groupslist)):
			subgroup(workpath, groupslist[i], sample_qiime, runsplit, "group")
			splitOTUtab(workpath, groupslist[i], runsplit)
			Vennplot(workpath, groupslist[i], subanalysis)
			alpha_plot(workpath, groupslist[i], subanalysis)
			pca_plot(workpath, groupslist[i], subanalysis)
			PCoA_NMDS_plot(workpath, groupslist[i], subanalysis)
			taxclass(workpath, groupslist[i], subanalysis, heatmap_width, heatmap_height)
			cluster_analysis(workpath, groupslist[i], subanalysis, cl_label_size, clbar_tree_height, blank_width, clbar_labelsize)
			LEfSe_DOTU(workpath, groupslist[i], subanalysis)
			myspeciestree(workpath, groupslist[i], subanalysis)
			PDFbuildgroup(workpath, groupslist[i])
			if len(environment) != 0:
				EnvironRDA(workpath, groupslist[i])
	if len(subgroups.split(";")[0]) != 0:
		subgroupslist = subgroups.split(";")
		for i in xrange(len(subgroupslist)):
			subgroup(workpath, subgroupslist[i], sample_qiime, runsplit, "subgroup")
			splitOTUtab(workpath, subgroupslist[i], runsplit)
			Vennplot(workpath, subgroupslist[i], subanalysis)
			alpha_plot(workpath, subgroupslist[i], subanalysis)
			pca_plot(workpath, subgroupslist[i], subanalysis)
			PCoA_NMDS_plot(workpath, subgroupslist[i], subanalysis)
			taxclass(workpath, subgroupslist[i], subanalysis, heatmap_width, heatmap_height)
			cluster_analysis(workpath, subgroupslist[i], subanalysis, cl_label_size, clbar_tree_height, blank_width, clbar_labelsize)
			LEfSe_DOTU(workpath, subgroupslist[i], subanalysis)
			myspeciestree(workpath, subgroupslist[i], subanalysis)
			PDFbuildgroup(workpath, subgroupslist[i])
