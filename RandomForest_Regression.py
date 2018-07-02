#!/usr/bin/env python

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.cross_validation import train_test_split
from sklearn.metrics import r2_score
import sys
import argparse
from argparse import RawTextHelpFormatter
import pickle
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt

def readparas(args):
	parser = argparse.ArgumentParser(description="The pipeline was designed to run the RandomForestRegression with sklearn.\n"
		"The developer: ShuangbinXu.\n"
		"Email:xusbin@anjiemed.com\n", formatter_class=RawTextHelpFormatter)
	parser.add_argument('Feature_tab',metavar="Feature_file", type=str, help="the Feature file, the format is:\n"
		"Feature\tF1\tF2\tF3\n"
		"sample1\t22\t3\t42\n")
	parser.add_argument('target_tab',metavar='target_file', type=str, help="the configure file, the format is:\n"
		"sample\tscore\n"
		"sample1\t89\n"
		"sample2\t76\n")
	parser.add_argument("-o", dest="outmodel", metavar="out_model", type=str, default="model_file.fit", help="the model file of RandomForestRegressor,\n"
		"default is model_file.fit\\")
	parser.add_argument("-T", dest="cpu", metavar="cpu_nums", type=int, default = -1, help="the cpu numbers of cross validation,\n"
		"default is -1 (all).")
	parser.add_argument("-x", dest="fold", metavar="fold_value", type=int, default = 10, help="the numbers of the fold in cross validation, default is 10")
	parser.add_argument("-p", dest="figure", metavar="out_figure", type=str, default = "figure_regression.svg", help="the figure about predict values and real values, default is figure_regression.svg")
	parser.add_argument("-t", dest="transpose", metavar="transpose", type=str, default="No", help="Transpose the Feature tab?, default is No")

        args = parser.parse_args()
        return vars(args)

def plotregression(y, predictresult, outfigure, R2, Spearman, Pearson):
	fig, ax = plt.subplots()
	ax.scatter(y, predictresult, edgecolors=(0, 0, 0))
	ax.plot([min(y), max(y)], [min(y), max(y)], "k--", lw=2)
	ax.annotate(R2, xy=(min(y), min(y)), xytext=(min(y), min(y)-1))
	ax.annotate(Spearman, xy=(min(y), min(y)), xytext=(min(y), min(y)-2))
	ax.annotate(Pearson, xy=(min(y), min(y)), xytext=(min(y), min(y)-3))
	ax.set_xlabel("Measured")
	ax.set_ylabel("Predicted")
	plt.savefig(outfigure, format="svg")


if __name__ == "__main__":
	params = readparas(sys.argv)
	feature = params["Feature_tab"]
	target = params["target_tab"]
	modelout = params["outmodel"]
	threas = params["cpu"]
	foldvalue = params["fold"] 
	figureout = params["figure"]
	Tflag = params["transpose"]
	feature = pd.read_table(feature, index_col=0, sep="\t", comment=None)
	if Tflag == "Yes":
		feature = feature.T
	target = pd.read_table(target, index_col=0, sep="\t", comment=None)
	target = target.fillna(target.median())
	tmptable = pd.merge(feature, target, left_index=True, right_index=True, how="outer")
	target = tmptable[target.columns.values].values.ravel()
	tmptable.drop(tmptable.columns[[-1]], axis=1, inplace=True)
        feature = tmptable

	#target = tmptable.values
	rfR = RandomForestRegressor(n_estimators=500, oob_score=True, random_state=100, n_jobs=threas)
	rfRscore = cross_val_score(rfR, feature, target, cv=foldvalue, n_jobs=threas)
	rfpredict = cross_val_predict(rfR, feature, target, cv=foldvalue, n_jobs=threas)
	rftest_R2score = r2_score(target, rfpredict)
	rfspearman = spearmanr(target, rfpredict)
	rfpearson = pearsonr(target, rfpredict)

	rfR.fit(feature, target)
	R2 = "Test data R-2 score: %s"%(round(rftest_R2score,3))
	Sp = "Test data Spearman correlation: %s"%(round(rfspearman[0],3))
	pe = "Test data Pearson correlation: %s"%(round(rfpearson[0],3))

	print ("Test data R-2 score: %s"%(round(rftest_R2score,3)))
	print ("Test data Spearman correlation: %s"%(round(rfspearman[0],3)))
	print ("Test data Pearson correlation: %s"%(round(rfpearson[0],3)))
	pickle.dump(rfR, open(modelout, "wb"))
	plotregression(target, rfpredict, figureout, R2, Sp, pe)
