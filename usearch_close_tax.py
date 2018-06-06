#!/usr/bin/env python
import sys

if len(sys.argv) != 4:
	print "#########################################################################"
	print "Function:The Script was designed to build the tax table of usearch close."
	print "Usage: python usearch_close_tax.py + taxonomy_txt + usearch_close_txt_tab + otuput"
	print "#########################################################################"
	sys.exit(1)
flagtmp = ["NR_074277.1", "JQ940228.1", "U55233.1", "X53217.1", "KJ567609.1.1488", "AF219235.1", "AF283539.2", "DQ831124.1.1414", "GU227149.1.1426", "KC771220.1.1511", "FR734082.1.1281", "CP014704.832081.833678", "HW313913.1.1460", "HW313943.1.1475", "HW313922.1.1489", "HW313923.1.1490", "HW313936.1.1440", "HW313951.1.1511", "HW313940.4.1583", "AZJT01000071.1165.2796", "AZTM01000094.703.2075", "AYSG01000002.292.2076", "EU547306.1.1464", "M58738.1", "AB269766.1", "CP012713.674346.675742", "AF070224.1", "EU118116.1", "FJ377885.1", "MH196340.1", "MG995565.1", "KF796667.1", "KY084558.1", "KC465400.1", "MF767595.1", "MF973087.1", "MG871232.1", "KX984112.1", "HQ658163.1", "EU734813.1", "KT964444.1", "NR_117620.1"]
fin = open(sys.argv[1], 'r')
tax = {}
for line in fin:
	tmp = line.replace("\n", "").split("\t")
	tax[tmp[0]] = tmp[1]
fin.close()

fin = open(sys.argv[2], "r")
fout = open(sys.argv[3], "w")
print >> fout,"# Constructed from biom file"
linenum = 0
for line in fin:
	tmp = line.replace("\n", "").split("\t")
	linenum += 1
	if linenum > 1:
		if tmp[0] in tax and tmp[0] not in flagtmp:
			print >> fout, "\t".join(tmp)+"\t"+tax[tmp[0]]
	else:
		tmp[0] = "#OTU ID"
		print >> fout, "\t".join(tmp)+"\ttaxonomy"
fin.close()
