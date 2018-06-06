#!/usr/bin/env python
import sys
import math

if len(sys.argv) != 2:
	print 'the program is designed to generate the annotation file for the graphlan to draw the species tree and distribution magic with the result of metaphlan.'
	print 'Usage:python generate_annotaion_for_species_tree.py + biom2lefseinputfile > output'
	sys.exit(1)

fin = open(sys.argv[1], 'r')
linenum = 0
sample = []
taxabun = {}
taxabun_o = {}
nump = 0
numo = 0
toptax = []
toptax_o = []
newtax = {}
for line in fin:
	linenum += 1
	tmp = line.replace('\n', '').split('\t')
	if linenum == 1:
		for i in xrange(1,len(tmp)):
			sample.append(tmp[i])
		newsample = list(set(sample))
	else:
		if tmp[0] in newtax:
			#newtax[tmp[0]]={}
			for i in xrange(len(sample)):
				if sample[i] in newtax[tmp[0]]:
					newtax[tmp[0]][sample[i]] += float(tmp[i+1])
				else:
					newtax[tmp[0]][sample[i]] = 0.0
					newtax[tmp[0]][sample[i]] += float(tmp[i+1])
		else:
			newtax[tmp[0]]={}
			for i in xrange(len(sample)):
				if sample[i] in newtax[tmp[0]]:
					newtax[tmp[0]][sample[i]] += float(tmp[i+1])
				else:
					newtax[tmp[0]][sample[i]] = 0.0
					newtax[tmp[0]][sample[i]] += float(tmp[i+1])

		taxtmp = tmp[0].split('|')
		if taxtmp[-1][0] == 'p':
			if tmp[0] in taxabun:
				for i in xrange(1,len(tmp)):
					taxabun[tmp[0]] += float(tmp[i])
			else:
				taxabun[tmp[0]]	= 0.0
				for i in xrange(1,len(tmp)):
					taxabun[tmp[0]] += float(tmp[i])
		elif taxtmp[-1][0] == 'f':
			if tmp[0] in taxabun_o:
				for i in xrange(1,len(tmp)):
					taxabun_o[tmp[0]] += float(tmp[i])
			else:
				taxabun_o[tmp[0]] = 0.0
				for i in xrange(1,len(tmp)):
					taxabun_o[tmp[0]] += float(tmp[i])
	taxpsort=sorted(taxabun.items(),key=lambda e:e[1],reverse=True)
	taxosort=sorted(taxabun_o.items(),key=lambda e:e[1],reverse=True)
fin.close()

for i in xrange(len(taxpsort)):
	nump += 1
	toptax.append(taxpsort[i][0])
	if nump >= 15:
		break
for i in xrange(len(taxosort)):
	numo += 1
	toptax_o.append(taxosort[i][0])
	if numo >= 15:
		break

fout = open(str(sys.argv[1])+'tmp2','w')
n = 0
for i in newtax:
	n += 1
	linetmp = []
	if n == 1:
		linetmp.append("ID")
	else:
		linetmp.append(i)
	for j in newtax[i]:
		if n == 1:
			linetmp.append(j)
		else:
			linetmp.append(str(newtax[i][j]))
	print >> fout, '\t'.join(linetmp)
fout.close()

annotlist = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
#samplecolors = ['#0000FF','#CC00FF','#FF0000','#FFA500','#800000','#483D8B','#006400','#8B2323','#33FF00','#A52A2A','#800080','#483D8B','#696969','#36648B']
#samplecolors = ["#C67CFF","#00BEC3","#7CAD00","#F8766E","#7CAD00", "#00BEC3", "#C67CFF"]
samplecolors = ['#00AED7','#FD9347', '#C1E168','#319F8C',"#F8766E","#7CAD00", "#00BEC3", "#C67CFF"]
taxcolors = ['#9ACD32','#EE6A50','#87CEFA','#FFC125','#D15FEE','#8DEEEE','#800000','#006400','#800080','#808080','#B0171F','#191970','#7B68EE','#FFC0CB','k']
tax2colors = {}
for i in xrange(len(toptax)):
	tax2colors[toptax[i]] = taxcolors[i]
sample2color = {}
print 'title	Species Tree And Distribution'
print 'title_font_size	33'
print 'total_plotted_degrees\t330'
print 'start_rotation\t270'
print 'branch_bracket_width\t0.0'
print 'annotation_background_alpha\t0.1'
print 'annotation_legend_font_size\t20'
print 'annotation_background_offset\t0.005'
print 'annotation_background_separation\t-0.03'
print 'annotation_background_width\t0.1'
print 'branch_thickness\t2.5'
print 'clade_marker_size\t20'
print 'clade_marker_edge_width\t0.8'
print 'class_legend_font_size\t25'
print 'class_legend_marker_size\t3'
for i in xrange(len(newsample)):
	print "ring_label_font_size\t"+str(i+1)+"\t25"
	print "ring_internal_separator_thickness\t"+str(i+1)+"\t0.7"
	print "ring_separator_color\t"+str(i+1)+'\t#888888'
	print "ring_label\t"+str(i+1)+"\t"+newsample[i]
	print "ring_label_color\t"+str(i+1)+"\t"+samplecolors[i]
	sample2color[newsample[i]] = samplecolors[i]
print "ring_label_font_size\t"+str(len(newsample)+1)+"\t20"
print 'ring_internal_separator_thickness\t'+str(len(newsample)+1)+'\t0.7'
print 'ring_separator_color\t'+str(len(newsample)+1)+'\t#888888'
print 'ring_label\t'+str(len(newsample)+1)+'\tAbundance when present'
print 'ring_label_color\t'+str(len(newsample)+1)+'\t#696969'


fin = open(str(sys.argv[1])+'tmp2', 'r')
linenum = 0
sample2=[]
for line in fin:
	linenum += 1
	tmp = line.replace('\n', '').split('\t')
	tmpFlag = line.replace('\n', '').split('\t')
	if linenum == 1:
		for i in xrange(1,len(tmp)):
			sample2.append(tmp[i])
	else:
	
		taxtmp = tmp[0].split('|')
		if taxtmp[-1][0] == "g":
			counttmp = 0.0
			tmpcoun = 0.0 
			del tmpFlag[0]
			tmpFlag = map(eval, tmpFlag)
			for i in xrange(1,len(tmp)):
				if float(tmp[i]) > counttmp:
					counttmp = float(tmp[i])
					j = i
				elif sum(tmpFlag) == 0.0:
					j = 1
				else:
					continue
			if j == "N" or j=="M":
				continue
			else:
				print tmp[0].split('|')[-1]+"\tring_color\t"+str(len(sample2)+1)+'\t'+sample2color[sample2[j-1]]
				print tmp[0].split('|')[-1]+"\tring_height\t"+str(len(sample2)+1)+'\t'+str(math.log(counttmp*5+0.002))
				print tmp[0].split('|')[-1]+"\tring_width\t"+str(len(sample2)+1)+"\t0.7"
			taxphly = taxtmp[0]+'|'+taxtmp[1]
			if taxphly in tax2colors:
				print tmp[0].split('|')[-1]+"\tclade_marker_color\t"+tax2colors[taxphly]
			else:
				continue
			for i in xrange(1,len(tmp)):
				tmpcoun += float(tmp[1])
			if float(tmpcoun) < 10.0:
				print tmp[0].split('|')[-1]+'\tclade_marker_size\t'+str(20)
			else:
				print tmp[0].split('|')[-1]+"\tclade_marker_size\t"+str(tmpcoun+20)
			print tmp[0].split('|')[-1]+'\tclade_marker_shape\t*'
			for i in xrange(len(sample2)):
				print tmp[0].split('|')[-1]+'\tring_color\t'+str(i+1)+"\t"+sample2color[sample2[i]]
				print tmp[0].split('|')[-1]+'\tring_height\t'+str(i+1)+'\t1.4'
				if float(tmp[i+1]) > 1.0:
					print tmp[0].split('|')[-1]+'\tring_alpha\t'+str(i+1)+'\t1.0'
				else:
					print tmp[0].split('|')[-1]+'\tring_alpha\t'+str(i+1)+'\t'+str(tmp[i+1])
fin.close()


	
for i in xrange(len(toptax_o)):
	print toptax_o[i].split('|')[-1]+"\tannotation\t"+annotlist[i]+": "+toptax_o[i].split('|')[-1]
	print toptax_o[i].split('|')[-1]+"\tannotation_background_color\tk"
	print toptax_o[i].split('|')[-1]+"\tannotation_font_size\t25"

for i in tax2colors:
	print i.split('|')[-1]+'\tclade_marker_color\t'+tax2colors[i]
	print i.split('|')[-1]+'\tclade_marker_size\t100'
	print i.split('|')[-1].split('__')[-1]+'\tclade_marker_color\t'+tax2colors[i]
	print i.split('|')[-1].split('__')[-1]+'\tclade_marker_size\t20'
