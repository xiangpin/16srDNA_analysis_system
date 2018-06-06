#!/usr/bin/env python

import sys
if len(sys.argv) !=3:
	print "#####################################"
	print "Function:The program is designed to extract the knowed_annotaion_tax tree"
	print "Usage:python pick_known_tax_to_tree.py metaphlan_out + g(p,c,o,f,g,s) > tmp.txt"
	print "#####################################"
	sys.exit(1)

fin = open(sys.argv[1], 'r')
a = sys.argv[2]
linenum = 0
groups = {}
tax = {}
classlist = ['k','p','c','o','f','g','s']
for line in fin:
	tmp = line.replace('\n', '').split('\t')
	linenum += 1
	if linenum > 1:
		taxtmp = tmp[0].split('|')
		p = classlist.index(a)	
		if taxtmp[-1][0] == "U" or  taxtmp[-1].split('_')[-1] == 'unclassified':
			continue
		elif classlist.index(taxtmp[-1][0]) <= p and len(taxtmp[-1]) > 3:
			taxname = '|'.join(taxtmp)
			tmp[0] = taxname
			linetmp = '\t'.join(tmp)
			print linetmp
		else:
			continue
			#linetmp = '\t'.join(tmp)
			#print linetmp
	else:
		linetmp = '\t'.join(tmp)
		print linetmp
		
fin.close()

