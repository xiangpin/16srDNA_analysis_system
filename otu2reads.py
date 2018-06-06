#!/usr/bin/env python

import sys

if len(sys.argv) !=2:
	print "Function: the program is designed to creat the otu2reads file with the usearch map.uc"
	print "Usage:python otu2reads.py results.uc"
	sys.exit(1)
uc = sys.argv[1]
otus = {}
with open(uc, "rt") as fin:
	for line in fin:
		tmp = line.replace("\n","").split("\t")
		flag = tmp[0]
		otu = tmp[9]
		seqid = tmp[8].split(";")[1].split("=")[1].strip()+"_"+tmp[8].split(";")[0].split(".")[len(tmp[8].split(";")[0].split("."))-1]
		if flag == "H":
			if otu not in otus:
				otus[otu] = []
				otus[otu].append(seqid)
			else:
				otus[otu].append(seqid)

fout = open("otu2seqs.tab","w")
for key in otus:
	tmp = ""
	for i in xrange(len(otus[key])):
		tmp += "\t"+otus[key][i]
	print >> fout , key  + tmp
fout.close()
