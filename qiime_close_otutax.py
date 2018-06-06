#!/usr/bin/env python
import sys
import os
import argparse


def mul_args(args):
	parser = argparse.ArgumentParser(description='the program is designed to change the tax of close qiime.')
	parser.add_argument('otu_tab',metavar='otu_tab',type=str,help='the otu_txt_tab of qiime output by close.')
	args = parser.parse_args()
	return vars(args)

def change(otu):
	classlist = ['k__','p__','c__','o__','f__','g__','s__']
	out="otu_tax_table.tab"
	fout = open(out,"w")
	with open(otu,'rt') as fin:
		linenum = 0
		for line in fin:
			linenum += 1
			tmp = line.replace("\n","").split("\t")
			if linenum >2:
				tax = tmp[-1]
				if tax == "Unassigned":
                                	#tmptax = "k__Unassigned; p__Unassigned; c__Unassigned; o__Unassigned; f__Unassigned; g__Unassigned; s__Unassigned"
                                	#tmp[-1] = tmptax
					continue
                        	else:
                                	taxlist = tax.split("; ")
                                	if len(taxlist) == 7:
                                        	for i in xrange(0,len(taxlist)):
                                                	if len(taxlist[i]) == 3:
                                                        	j = i	
                                                       		for ind in xrange(j, 7):
                                                                	taxlist[ind] = taxlist[ind]+"un_"+taxlist[j-1][0]+"_"+taxlist[j-1][3:]
                                	else:
                                        	for i in xrange(0, len(taxlist)):
                                                	if len(taxlist[i]) == 3:
                                                        	j = i
                                                        	for ind in xrange(j, len(taxlist)):
                                                                	taxlist[ind] = taxlist[ind]+"un_"+taxlist[j-1][0]+"_"+taxlist[j-1][3:]
                                        	for e in xrange(len(taxlist), 7):
                                                	n = len(taxlist) - 1
                                                	taxlist.append(classlist[e] + "un_" + taxlist[n][0]+"_"+taxlist[n][3:])
                                	tmptax = "; ".join(taxlist)
                                	tmp[-1] = tmptax
                        	tmpline = "\t".join(tmp)
                        	print >> fout, tmpline
			else:
				print >> fout, "\t".join(tmp)


if __name__ == "__main__":
	params = mul_args(sys.argv)
	otu_tab = params['otu_tab']
	change(otu_tab)					
