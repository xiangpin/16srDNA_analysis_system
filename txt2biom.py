#!/usr/bin/env python
import sys,os

#def myparse():
if len(sys.argv) !=3:
	print "Function:This program was designed to covert the otu_table_txt file to the otu_table_biom file."
	print "Usage:python txt2biom.py otu_tax_table.txt otu_tax_table.biom"
	sys.exit(1)
otu_txt_table = sys.argv[1]
otu_biom_table = sys.argv[2]
cmd = 'biom convert -i ' + otu_txt_table +' -o ' + otu_biom_table +' --table-type="Taxon table" --to-hdf5 --process-obs-metadata taxonomy'
os.system(cmd)
