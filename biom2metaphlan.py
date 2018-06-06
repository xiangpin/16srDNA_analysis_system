#!/usr/bin/env python
import Bio
import argparse as ap
import sys
import json
import csv
import numpy as np

reload(sys)
sys.setdefaultencoding('utf-8')
try:
    from biom.parse import parse_biom_table
    import biom.parse
    import biom.table
except:
    sys.stdout.write("Error, the BIOM library (http://biom-format.org/) is not installed\n")
    sys.exit()

__author__ = 'Sami Pietila (sampie@iki.fi)'
__version__ = '0.8'
__date__ = '27th Sep 2013'


class TaxonNode:
    def __init__(self): 
        self.parent = None
        self.fracs = None # one for each sample
        self.counts = None # one for each sample
        self.taxon = None
        self.children = None
        self.include = False
        self.sampleIDs = None # defined only for root
        self.maxdepth = 0 # defined only for root
        
#     def hasNode(self, taxon):
#         if self.children != None:
#             for child in self.children:
#                 if child.taxon == taxon:
#                     return True
#         return False
        
    def addNode(self, taxon, fracs, counts):
        tn = TaxonNode()
        tn.parent = self
        tn.taxon = taxon
        tn.fracs = fracs
        tn.counts = counts
        if (self.children == None): 
            self.children = []
        self.children.append(tn)
        return tn
        
    def getNode(self, taxon):
        if self.children != None:
            for child in self.children:
                if child.taxon == taxon:
                    return child
        return None    
         
                     
    def printTree(self):
        zList = [self]
        while len(zList)>0:
            newL = []
            for l in zList:
                if l.children == None:
                    c = l
                    while (c != None):
                        if c.taxon != None:
                            sys.stdout.write('|')
                            sys.stdout.write(c.taxon)
                            sys.stdout.write("(")
                            if c.fracs != None:
                                for fr in c.fracs:
                                    sys.stdout.write(str(fr))
                                    sys.stdout.write(",")
                            if c.counts != None:
                                for co in c.counts:
                                    sys.stdout.write(str(co))
                                    sys.stdout.write(",")
                            sys.stdout.write(")")
                        c = c.parent
                    sys.stdout.write('\n')
                
                if l.children != None:
                    for child in l.children:
                        newL.append(child)
            zList = newL


    # run for root node        
    def countsToUnclassifiedChild(self):
        taxonIDs = ["k__","p__","c__","o__","f__","g__", "s__"]

        zList = [self]
        while len(zList)>0:
            newL = []
            for l in zList:
                if l.children != None:
                    for child in l.children:
                        newL.append(child)
                    if l.counts != None:
                        if (l.taxon != None):
                            #print l.taxon, ", l.depth():", l.depth()
                            postfix = "_unclassified"
                            i = taxonIDs.index(l.taxon[0:3])
                            taxon = taxonIDs[i+1]+l.taxon[3:]+postfix
                        else:
                            taxon = "k__unclassified"
                        l.addNode(taxon, None, l.counts)
                        l.counts = None

            zList = newL
            fracs = None

    # run for root
    def nodesByLevel(self):
        zList = [self]
        zLists = []
        while len(zList)>0:
            zLists.append(zList)
            newL = []
            for l in zList:
                if l.children != None:
                    newL.extend(l.children)
            zList = newL
        return zLists

    
    def getLeafNodes(self):
        leafNodes = []
        zList = self.children
        while len(zList)>0:
            newL = []
            for l in zList:
                if l.children == None:
                    leafNodes.append(l)
                else:
                    newL.extend(l.children)
            zList = newL
        return leafNodes

    def depth(self):
        #taxonIDs = ["k__","p__","c__","o__","f__","g__", "s__"]
        c = self.parent
        depth = 0
        while c != None and c.parent != None:
            c = c.parent
            depth += 1
        return depth    

    # run to root. if leaf nodes are not at max level, extend them.
    def extendLeafNodes(self):
        taxonIDs = ["k__","p__","c__","o__","f__","g__", "s__"]
        leafNodes = self.getLeafNodes()
        for l in leafNodes:
            c = l
            if c.taxon[0:3] in taxonIDs:
                taxonbody = c.taxon[3:]
            else:
                taxonbody = c.taxon
                
            if taxonbody[-(len("_unclassified")):] == "_unclassified":
                taxonbody = taxonbody[0:-len("_unclassified")]
            depth = l.depth()
            while depth < root.maxdepth:
                postfix = "_unclassified"
                taxon = taxonIDs[depth+1]+taxonbody+postfix
                newNode = c.addNode(taxon, None, None)
                newNode.counts = c.counts
                c.counts = None
                depth += 1
                c = newNode
            #print l.taxon,":", l.depth()
            

    # run for root node        
    def propagateFracs(self):
        zLists = self.nodesByLevel()
        zLists.reverse()
        for zList in zLists:
            zTotal = [0]*len(self.sampleIDs)
            for node in zList:
                if node.parent != None:
                    if node.parent.counts == None:
                        node.parent.counts = node.counts
                    else:
                        node.parent.counts = map(sum, zip(node.parent.counts,node.counts))
                zTotal = map(sum, zip(zTotal,node.counts))
#            print "zTotal:", zTotal[0]
            for node in zList:
                node.fracs = [float(0)]*len(self.sampleIDs)
                for i in range(len(zTotal)):
                    if zTotal[i] > 0:
                        node.fracs[i] = float(node.counts[i]) / float(zTotal[i])
#                        if j == 7 and i == 0:
#                            print node.counts[i], ":", zTotal[i], ":", float(node.counts[i]) / float(zTotal[i])

        
    def taxaList(self):
        taxaL = []
        n = self
        while n != None and n.parent != None: # Root does not have taxa
             taxaL.append(n.taxon)
             n = n.parent
        taxaL.reverse()     
        return taxaL


def readBIOM(fileName):
    f = open(fileName, "r")
    table = parse_biom_table(f)
    #print table
    f.close()

    root = TaxonNode()

    root.sampleIDs = list(table._sample_ids)
    #print "SampleIDs type:", type(table.SampleIds)

    #for obs in table.iterObservations():
    for obs in table.iter(axis='observation'):
	
        counts = obs[0]
        otuName = obs[1]
        taxonomy = obs[2]["taxonomy"]
        root.maxdepth = max(root.maxdepth,len(taxonomy)-1)
        # Build Tree
        node = root
        for taxon in taxonomy:
            n = node.getNode(taxon)
            if n == None:
                node = node.addNode(taxon, None, None)
            else:
                node = n
            if taxon == taxonomy[-1]:
                if node.counts != None:
                    node.counts = map(sum, zip(node.counts, counts))
                else:
                    node.counts = counts
    return root


def taxaAggregate(taxaList, targetLevel):
    level = len(taxaList)-1
    taxonL = []
    while level > targetLevel:
        taxonL.append(taxaList[level])
        level -= 1
    #taxonL.reverse()
    if len(taxonL) > 0:
        return "/".join(taxonL)
    else:
        return taxaList[-1]


# Write samples in a single file
class singleWriter:
    def __init__(self, sampleHeaders, taxaLevel):
        self.file = None
        self.sampleHeaders = sampleHeaders
        self.taxaLevel = taxaLevel
    def open(self, filename, root):
        self.file = open(filename, "w")
        self.file.write("\t".join(["ID"]+self.sampleHeaders) + "\n")
    def writeSample(self, node):
        taxaL = node.taxaList()
        if self.taxaLevel != None:
            aggregatedTaxa = taxaAggregate(taxaL, self.taxaLevel)
            if aggregatedTaxa != None:
                taxaL[-1] = aggregatedTaxa
        taxonomyString = "|".join(taxaL)
        fracsString = "\t".join(str(f*100) for f in node.fracs)
        self.file.write(taxonomyString + "\t" + fracsString + "\n")
    def close(self):
        self.file.close()
        

# Write samples in separate files
class multiWriter:
    def __init__(self, taxaLevel):
        self.files = {}
        self.root = None
        self.taxaLevel = taxaLevel
    def open(self, filenamepostfix, root):
        self.root = root
        for sID in root.sampleIDs:
            f = open(sID + filenamepostfix, "w")
            self.files[sID]=f
#            f.write("ID\t" + sID + "\n")
    def writeSample(self, node):
        for i, sID in enumerate(self.root.sampleIDs):
            taxaL = node.taxaList()
            if self.taxaLevel != None:
                aggregatedTaxa = taxaAggregate(taxaL, self.taxaLevel)
                if aggregatedTaxa != None:
                    taxaL[-1] = aggregatedTaxa
            taxonomyString = "|".join(taxaL)
            fracsString = str(node.fracs[i]*100)
            self.files[sID].write(taxonomyString + "\t" + fracsString + "\n")
    def close(self):
        for sID in self.root.sampleIDs:
            self.files[sID].close()


def writeMetaPhLan(writer, filename, root):
    writer.open(filename, root)
    node = root.children[0]
    index = 0
    while node != None:
        if index == 0:
            writer.writeSample(node)

        if node.children != None:
            l = len(node.children)
        
        if node.children != None and index < len(node.children):
            node = node.children[index]
            index = 0
        else:
            if node.parent != None:
                index = node.parent.children.index(node)+1
            node = node.parent

    writer.close()
    return


def read_params(args):
    p = ap.ArgumentParser( description= 
                           
            "DESCRIPTION\n"
            "  metaphlan2biom.py version "+__version__+" ("+__date__+")\n" 
            
            "  metaphlan2biom.py is a tool to convert Qiime produced\n"
            "  OTU table to MetahPhlAn's relative abundance profile of metagenome.\n" 
            "  AUTHOR: "+__author__+"\n\n"
            "EXAMPLE\n"
            "  python biom2metaphlan.py -i otu_table.biom -o otu_table_metahplan.txt\n"
            "  for write samples in separate files (for krona):\n"
            "  python biom2metaphlan.py -i otu_table.biom -o .txt -s\n"
            "  for write lefse output using metadata from qiime mapping file:\n"
            "  python biom2metaphlan.py -l -m map-file.txt -g group-name -i otu_table.biom -o metahplan.txt\n",
            formatter_class=ap.RawTextHelpFormatter )
    arg = p.add_argument
    arg( '-i', metavar='INPUT_FILE', type=str, nargs='?', default=None, help= 
         "Input file (OTU table in BIOM format).\n")
        
    arg( '-o', metavar='OUTPUT_FILE', type=str, nargs='?', default=None, help= 
         "Output file to be written. The output file in MetaPhlAn txt format.\n")

    arg( '-m', metavar='METADATA_FILE', type=str, nargs='?', default=None, help= 
         "Qiime mapping file to read metadata.\n")

    arg( '-g', metavar='GROUP', type=str, nargs='?', default=None, help= 
         "Group ID from Qiime metadata file to use to print Lefse Output.\n")

    arg( '-a', metavar='TAXA_LEVEL', type=int, nargs='?', default=None, help= 
         "Aggregate taxa annotation to given level (4=genus, 3=family, 2=order, 1=class, 0=phylum).\n")

    arg( '-l', action='store_true', help= 
         "Write Lefse format. Also -m and -g needed.\n")

    arg( '-s', action='store_true', help= 
         "Samples are being written in separate files. \n"
         "(-o gives postfix for sample files.) \n")
            
    return vars(p.parse_args())


def readQiimeMetaDataFile(fileName):
    f = open(fileName, "r")
    fileData = f.read()
    f.close()
    
    lines = fileData.split("\n")
    header = lines[0].split("\t")
    header[0] = header[0][1:]
    headerLen = len(header)
    
    mdata = [header]
    
    for j in range(len(lines)):
        if lines[j] == "" or lines[j].startswith("#"): 
            continue
        line = lines[j].split("\t")
        if len(line) != headerLen: 
            raise Exception("A data line parameter count does not match header: ", lines[j])
        mdata.append(line)

    return mdata


# Generate list of groups for samples
def groupList(metadata, groupName, sampleIDs):
    sampleGroups = {}
    header = metadata[0]
    
    for l in metadata[1:]:
        sampleGroups[l[0]] = l[header.index(groupName)]
    
    groupList = []
    for sample in sampleIDs:
        groupList.append(sampleGroups[sample])
        
    return groupList            
        
if __name__=="__main__":
    pars = read_params( sys.argv )
    metaData = None

    if pars['i'] is None:
        raise Exception("Use -i to give OTU table in BIOM format.")
    if pars['o'] is None:
        raise Exception("Use -o to give file name for MetaPhlAn txt format file.")
    if pars['l']:
        if pars['m'] is None:
            raise Exception("Use -m <mapping-file-name> to qime Qiime metadata (mapping) file.")
        if pars['g'] is None:
            raise Exception("Use -g <group name> to give group name (in metadata file) by which identify biomarkers.")
        metaData = readQiimeMetaDataFile(pars['m'])
        if pars['g'] not in metaData[0]:
            raise Exception(groupName," group is not found from metadata file.")
 
    root = readBIOM(pars['i'])
    root.countsToUnclassifiedChild()
    root.extendLeafNodes()
    root.propagateFracs()

    writer = None 
    if pars['l']:
        groupList = groupList(metaData, pars['g'], root.sampleIDs)
        writer = singleWriter(groupList, pars['a'])
    else:
        if pars['s']:
            writer = multiWriter(pars['a'])
        else:
            writer = singleWriter(root.sampleIDs, pars['a'])
    
    writeMetaPhLan(writer, pars['o'], root)
        
    
    
    
 
