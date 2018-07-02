#!/usr/bin/env python

import sys
import argparse
from argparse import RawTextHelpFormatter
import pandas as pd
import numpy as np
from scipy import interp
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_val_score
from sklearn import preprocessing
from sklearn.preprocessing import LabelBinarizer
from sklearn.model_selection import StratifiedKFold

def readparas(args):
	parser = argparse.ArgumentParser(description="The Script was designed to run the RandomForestClassifier with sklearn.\n"
	"The developer: ShuangbinXu.\n"
	"Email:xusbin@anjiemed.com\n", formatter_class=RawTextHelpFormatter)
	parser.add_argument('Feature_tab',metavar="Feature_file", type=str, help="the Feature file, the format is:\n"
		"header line or No\n"
		"sample\tsample1\tsample2\tsample3\tsample4\n"
		"Feature1\t34\t242\t2\t23\n")		
	parser.add_argument('group_tab',metavar='group_file', type=str, help="the configure file, the format is:\n"
		"sample\tgroup\n"
		"sample1\tclass1\n"
		"sample2\tclass2\n"
		"sample3\tclass3\n")
	parser.add_argument("-o", dest="outmodel", metavar="out_model", type=str, default="model_file.fit", help="the model file of RandomForestClassifier,default is model_file.fit")
	parser.add_argument("-T", dest="cpu", metavar="cpu_nums", type=int, default = -1, help="the cpu numbers of cross validation, default is -1 (all).")
	parser.add_argument("-x", dest="fold", metavar="fold_value", type=int, default = 10, help="the numbers of the fold in cross validation, default is 10")
	#parser.add_argument("-t", dest="transpose", metavar="transpose", type=str, default="No", help="Transpose the Feature tab?, default is No")
	parser.add_argument("-p", dest="figure", metavar="out_figure", type=str, default = "figure_AUC_ROC.svg", help="the ROC figure, default is figure_AUC_ROC.svg")
	args = parser.parse_args()
	return vars(args)


if __name__ == "__main__":
	params = readparas(sys.argv)
	otutab = params["Feature_tab"]
	samplefile = params['group_tab']
	threas = params["cpu"]
	foldvalue = params["fold"]
	#Tflag = params["transpose"]
	modelout = params["outmodel"]
	figureout = params["figure"]	

	otu = pd.read_table(otutab, skiprows=0, index_col=0, header=1, sep="\t", comment=None)
	otu.pop("taxonomy")
	otu = otu.T
	sample = pd.read_table(samplefile, index_col=0, header=0, sep="\t", comment=None)
	tmptable = pd.merge(sample, otu, left_index=True, right_index=True, how="outer")
	target = np.array(tmptable.pop("group"))
	label = np.unique(target)
	otu = pd.DataFrame(tmptable).values
	le = preprocessing.LabelBinarizer()
	
	target = le.fit_transform(target)
	nclass = target.shape[1]
	cv = StratifiedKFold(n_splits=foldvalue)
	colorlist = ["#C1E168", "#319F8C", "#FD9347", "#00AED7"]
	rf = RandomForestClassifier(n_estimators = 1000, oob_score=True, n_jobs=threas, random_state=1000)
	plt.figure(1,figsize=(5,4.9))
	for i in range(nclass):
		tmplabel = label[i]
		tmpcolor = colorlist[i] 
		y = target[:,i]
		tprs = []
		aucs = []
		mean_fpr = np.linspace(0, 1, 100)
		for train, test in cv.split(otu, y):
			probas_ = rf.fit(otu[train], y[train]).predict_proba(otu[test])
			fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
			tprs.append(interp(mean_fpr, fpr, tpr))
			tprs[-1][0] = 0.0
			roc_auc = auc(fpr, tpr)
			aucs.append(roc_auc)
		mean_tpr = np.mean(tprs, axis=0)
		mean_tpr[-1] = 1.0
		mean_auc = auc(mean_fpr, mean_tpr)
		std_auc = np.std(aucs)
		plt.plot(mean_fpr, mean_tpr, lw=1.8, alpha=0.8, color=tmpcolor,
			label = 'Class {0} Mean ROC (AUC = {1:0.2f} $\pm$ {2:0.2f})'
			''.format(tmplabel, mean_auc, std_auc)
			)
	plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='black',
	         label='Luck', alpha=.6)
	plt.xlabel('False positive rate', fontsize=10)
	plt.ylabel('True positive rate',fontsize=10)
	plt.xticks(size = 7)
	plt.yticks(size = 7)
	plt.title('ROC curve of the class\n', fontsize=13, fontweight="bold")
	plt.legend(loc='lower right', fontsize=8)
	
	plt.savefig(figureout, format="svg")
