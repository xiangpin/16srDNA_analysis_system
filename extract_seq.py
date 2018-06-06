#!/usr/bin/env python
import sys

if len(sys.argv) !=3:
	print "###################################################################################"
	print "Function:The script was designed to extract the repsent sequence from the database."
	print "Usage:python extrat_seq.py + database + otu_tab."
	print "###################################################################################"
	sys.exit(1)

fin = open(sys.argv[1], 'r')
seq = {}
for line in fin:
	tmp = line.replace("\n", "")
	if tmp.startswith(">"):
		key = tmp
		seq[key] = ""
	else:
		seq[key] += tmp
fin.close()

fin = open(sys.argv[2], 'r')
for line in fin:
	tmp = line.replace("\n", "").split("\t")
	if ">"+tmp[0] in seq:
		print ">"+tmp[0]+'\n'+seq[">"+tmp[0]]
fin.close()
