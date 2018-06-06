#!/usr/bin/env python
import sys

if len(sys.argv) != 4:
	print "#######################################"
	print "Function:The script was used to build the function plot input"
	print "Usage:python PICRUSt_out_change.py + cog_predict.txt + output +type_function(cog,ko)"
	print "#######################################"
	sys.exit(1)
fin = open(sys.argv[1], 'r')
fout = open(sys.argv[2], 'w')
type_function = sys.argv[3]
linenum = 0
for line in fin:
	linenum += 1
	linetmp = []
	tmp = line.replace('\n', '').split('\t')
	if linenum == 2:
		linetmp.append('OTU ID')
		for i in xrange(1,len(tmp)-1):
			linetmp.append(tmp[i])
		print >> fout, '\t'.join(linetmp)
	elif linenum > 2:
		if type_function == "cog":
			linetmp.append(tmp[0]+"@"+tmp[-1])
			for i in xrange(1,len(tmp)-1):
				linetmp.append(tmp[i])
			print >> fout, '\t'.join(linetmp)
		else:
			for i in xrange(0,len(tmp)-1):
				linetmp.append(tmp[i])
			print >> fout, '\t'.join(linetmp)
fin.close()
fout.close()
