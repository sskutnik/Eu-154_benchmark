import matplotlib 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import re
import itertools

sns.set(style="ticks",font_scale=1.5)
sns.set_context(rc={'lines.markeredgewidth': 1.0, 'lines.markeredgecolor': 'k', 'lines.linecolor': 'k'})

IMAGE_DEST='./images/'
TEX_DEST='./latex/'
GLOB_COL_WRAP = 2 # number of columns per line for small multiples
OUTPUT_FORMAT = 'pdf'
SHOW_PLOTS = False
DEFAULT_YLIM = [0.4,1.3] # General y-axis limits for all plots

# Isotope regular expression, e.g. "233U"
RX_ISO = re.compile(r"([0-9]+)([a-zA-Z]+)")
# Alternate isotope regex, e.g. "U233" / "U-233"
RX_ISO2 = re.compile(r"([a-zA-Z]+)[-]?([0-9]+)")
# Properly-formatted element (U-233)
RX_ISO3 = re.compile(r"([a-zA-Z]+)-([0-9]+)")

# Calculate C/E from a calculations dataframe composed of multiple evaluations
# and a measurements dataframe, then add C/E back to the calculations dataframe
def calculateCE(df_calc, df_measured, evalUnits):

	CE_series = []
	sigma_series = []
	# Match measurements to calculations by sample
	for sample in df_calc["Sample"].unique():
		df_measSample = df_measured[df_measured["Sample"] == sample].reset_index()
		
		# Loop over individual evaluations within samples for C/E values
		for eval in df_calc["Evaluation"].unique():
			isEval = df_calc[(df_calc["Evaluation"] == eval) & (df_calc["Sample"] == sample)].reset_index()
			tmpSeries = pd.Series(isEval[evalUnits] / df_measSample["Inventory"], name = "C/E")
			CE_series.append(tmpSeries)
			sigma_series.append(pd.Series(df_measSample["Uncertainty"], name = "Sigma"))
			
	# Add the C/E & measurement uncertainty data back to the measurements			
	joinedCE = pd.concat(CE_series, ignore_index=True)
	joinedSigma = pd.concat(sigma_series, ignore_index=True)
	df_calc["C/E"] = joinedCE
	df_calc["Sigma"] = joinedSigma 
	
	#df_calc = df_calc.join(joinedCE)
	#print(df_calc.keys())
	
# Clean up isotopes - e.g., from "233U" and "233u" to "U-233"
def cleanupIsoNames(df_iso):
	isoNames = []
	eleNames = []
	for line in df_iso["Isotope"]:
		element = ""
		if(RX_ISO.match(line)):
			element = re.sub(RX_ISO,r'\2',line).title()
			line = re.sub(RX_ISO,r'\2-\1',line).title()
		elif(RX_ISO2.match(line)):
			element = re.sub(RX_ISO2,r'\1',line).title()
			line = re.sub(RX_ISO2,r'\1-\2',line).title()
	
		eleNames.append(element)
		isoNames.append(line)

	df_iso["Isotope"] = pd.Series(isoNames)
	df_iso["Element"] = pd.Series(eleNames)

	
def plotIsos(df_iso, labels, seriesInfo, pltTitle, filename, sigmaLabel="Sigma", yLim=DEFAULT_YLIM, legendLoc="best",burnup=0.0):
	
    for idx in range(0,len(labels)):
	    # Using .values attribute to handle cases where we exclude the first element (e.g., U-233), which plt.errorbar tries to access
        plt.errorbar(df_iso["Isotope"].index, df_iso[labels[idx]].values,df_iso[sigmaLabel].values,ls='',**seriesInfo[idx])

    if(burnup > 0):
        pltTitle = pltTitle = "\nBurnup: {:.2f} GWd/MTU"
		
    plt.xticks(df_iso.index,df_iso["Isotope"])
    plt.xlim(min(df_iso.index)-1,max(df_iso.index)+1)
    plt.title(pltTitle)
    plt.ylabel("C/E")
    plt.ylim(yLim)
    #plt.ylim(bottom=0.8)
    sns.despine(top=True)
    plt.legend(loc=legendLoc)
    plt.tight_layout()
    if SHOW_PLOTS:
        plt.show()
    else:
        plt.savefig(IMAGE_DEST + filename + '.' + OUTPUT_FORMAT)
        plt.close()

def plotIsosFP(df_iso, assemName, seriesInfo, filename, burnups={},plotAspect=1.05,colWrap=GLOB_COL_WRAP, fontScale=1.5):
	
	sns.set(style="ticks",font_scale=fontScale)
	
	isoNumbers = []
	idx = 0
	for isotope in df_iso["Isotope"].unique():
		isoNumbers.append((isotope, idx))
		idx = idx + 1
	isoDict = dict(isoNumbers)
	
	isoMarkers = ['o', 'o', 's', '^', 'v', '*', 'H']
	edgeColors = ['k'] * len(df_iso['Evaluation'].unique())
	edgeWidths = ['1'] * len(df_iso['Evaluation'].unique())
	df_iso["IsoNum"] = df_iso["Isotope"].map(isoDict)
	
	palette = itertools.cycle(sns.color_palette())
	ebDict = { }
	
	for eval in df_iso['Evaluation'].unique():
		ebDict[eval] = palette.next()

	g = sns.factorplot(data=df_iso,kind='point',col='Sample',x='Isotope',y='C/E',hue='Evaluation',legend=False,scale=1.75,size=4,
	                   col_wrap=colWrap,linestyles='None',aspect=plotAspect,markers=isoMarkers)

	for ax in g.axes:
		for mk in ax.collections:
			mk.set_edgecolor('black')
			mk.set_linewidth(0.075)

	g.map_dataframe(facetErrorBars,"IsoNum","C/E","Sigma",ebDict,data=df_iso, hue='Evaluation',zorder=0,elinewidth=1.5,markeredgecolor='k',ls='')

	g.set_xlabels('Isotope').set_titles("{col_name}")
	plt.legend(bbox_to_anchor=(1.1, 0.90), loc=2, borderaxespad=0., title="Evaluation")
	sns.despine(top=True)
	
	titles = []
	for sample in df_iso["Sample"].unique():
		tmpTitle = "Sample {:s}".format(sample)
		if(sample in burnups):
			tmpTitle = tmpTitle + "\nBurnup: {:.2f} GWd/MTU".format(burnups[sample])
		titles.append(tmpTitle)

	for ax, plotTitle in zip(g.axes.flat, titles):
		ax.set_title(plotTitle)
	
	
#	plt.legend(loc=legendLoc)
	plt.tight_layout()
	if SHOW_PLOTS:
		plt.show()
	else:
		plt.savefig(IMAGE_DEST + filename + '.' + OUTPUT_FORMAT)
		plt.close()

def facetErrorBars(x, y, sigma, cDict, color=None, label=None, **kwargs):
	# Add error bars to a FacetGrid / FactorPlot and match the bar color to the hue
	data = kwargs.pop("data")
	hue = kwargs.pop("hue")

	ebColors = []
	for eval in data[hue].unique():
		ebColor = cDict[eval]
		xVals = data[x].loc[data[hue] == eval].values
		yVals = data[y].loc[data[hue] == eval].values
		yErrs = data[sigma].loc[data[hue] == eval].values	
		plt.errorbar(xVals, yVals, xerr=None, yerr=yErrs, ecolor=ebColor, color=None,**kwargs)
		
def plotIsoBuFP(df_iso, labels, seriesInfo, isoName, sigmaLabel="Sigma", plotAspect=1.05, plotFit=False):

	plt.clf()
	isoMarkers = ['o', 'o', 's', '^', 'v', '*', 'H']

	for idx in range(0,len(labels)):
	#	print("Label = ", labels[idx])
		# Need to down-select either sample or evaluation as a series - plot C/E values for this...
		df_plot = df_iso.loc[(df_iso["Evaluation"] == labels[idx])]
		eBar = plt.errorbar(df_plot["Burnup"].values, df_plot["C/E"].values,df_plot[sigmaLabel].values,marker=isoMarkers[idx],label=labels[idx],ls='',markeredgecolor='k', markeredgewidth=0.15, ms=9)
		if(plotFit):
			eColor = eBar[0].get_color()
			fitParams = np.polyfit(df_plot["Burnup"].values, df_plot["C/E"].values, 1, w=1.0/df_plot["Sigma"].values, full=False)

			
			buVals = df_plot["Burnup"].values
			ybar = np.mean(df_plot["C/E"].values)
			ssreg = np.sum(((np.poly1d(fitParams)(buVals)-ybar)/df_plot["Sigma"])**2)
			sstot = np.sum(((df_plot["C/E"].values-ybar)/df_plot["Sigma"])**2)
			#print(ssreg,sstot)
			r2 = ssreg/sstot
			print("Evaluation: ", labels[idx], "Slope: {:.4f}, {:.4f}, R^2={:.4f}".format(fitParams[0], fitParams[1], r2))

			buVals.sort()
			plt.plot(buVals, np.poly1d(fitParams)(buVals), color=eColor,lw=0.75,ls='--')
			
	sns.despine(top=True)
	plt.xlabel("Burnup (GWd/MTU)")
	plt.ylabel("C/E")
	plt.title(isoName + " C/E variation with burnup")
	plt.xlim(np.floor(np.min(df_iso["Burnup"].values)-2.0), np.ceil(np.max(df_iso["Burnup"].values))+2.0)
	
	lgd = plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left', borderaxespad=0., title="Evaluation")

	filename = isoName + "_CE_BU"
	
	#
	if SHOW_PLOTS:
		plt.tight_layout()
		plt.show()
	else:
		#plt.savefig(IMAGE_DEST + filename + '.' + OUTPUT_FORMAT)
		plt.savefig(IMAGE_DEST + filename + '.' + OUTPUT_FORMAT,bbox_extra_artists=(lgd,), bbox_inches='tight')
		plt.close()

	
def generateLatexTables(df_iso, df_measured, filename):
	#texTable = "\\begin{table}[htb] \n\\centering \n"
	for sample in df_iso["Sample"].unique():
		evalCount = len(df_iso["Evaluation"].unique())
		texTable = "% \\begin{longtable}{l l | " + "m{1.3cm} " * 2 + "m{2.8cm} " * (evalCount - 2) + "}\n"
		texTable += "%\n% CAPTION GOES HERE \n%\n"
		texTable += "\\multirow{2}{*}{Isotope} & \\multirow{2}{*}{$\sigma_{meas.}$} & "
		texTable += "\\multicolumn{{{:d}}}{{c}}{{C/E}} \\\\ \n & ".format(evalCount)	
		
		for eval in df_iso["Evaluation"].unique():
			texTable += " & " + eval
		
		texTable += "\\\\ \n \\toprule \n"
		texTable += "\\endhead \n"
		
		# Write out C/E for each isotope (rows) / evaluation (columns)
		for isotope in df_iso["Isotope"].unique():
			isoEvals = df_iso.loc[(df_iso["Isotope"] == isotope) & (df_iso["Sample"] == sample)]
			measUncert = df_measured.loc[(df_measured["Isotope"] == isotope) & (df_measured["Sample"] == sample)]["Uncertainty"]
	
			# Only print isotope if any of the evaluations is not NaN
			if(isoEvals["C/E"].notnull().values.any()):
				sigmaString = ""
				if(measUncert.notnull().all()):
					sigmaString = " & {:.4f} ".format(measUncert.item()) 
				else:
					sigmaString = " & --- "
					
				texTable += isoToLatex(isotope) + sigmaString
				for index, row in isoEvals.iterrows():
					#print(row['C/E'])
					texTable += "& {:.4f} ".format(row['C/E'])
	
				texTable += "\\\\ \n"
		
		texTable += "\\bottomrule \n"
		texTable += "%\\end{longtable}\n"
		texFile = open(TEX_DEST + filename + "_" + sample + ".tex", "w")
		texFile.write(texTable)
		texFile.close()
		#print("SAMPLE: ", sample, "\n", texTable)
			#tmpTitle = "Sample {:s}".format(sample)

# Convert a properly-formatted isotope (U-233) to LaTeX syntax (${}^{233}{U}$)
def isoToLatex(isoName):
	latexIso = ""
	if (RX_ISO3.match(isoName)):
		latexIso = RX_ISO3.sub(r'${}^{\2}\mathrm{\1}$',isoName)
	#print(latexIso)
	return latexIso


#================================================
# Measurement data 
#===============================================
	
labels = ["ENDF VII.0","ENDF VII.1","ENDF VII.0 + mod-Wright","ENDF VII.0 + mod-Mughabghab", "ENDF VII.1 + mod-Mughabghab", "JEFF-3.2 + mod-Mughabghab"]
# Default series meta-info
defaultSeriesInfo = [{ 'markersize':8,'marker':'o','label':'ENDF VII.0 (nominal)' },
                     { 'markersize':8,'marker':'o','label':'ENDF/VII.1 (nominal)' },
                     { 'markersize':8,'marker':'^','label':'VII.0 + mod. Wright' },
                     { 'markersize':8,'marker':'v','label':'VII.0 + mod. Mughabghab' },
                     { 'markersize':8,'marker':'>','label':'VII.1 + mod. Mughabghab' },
					 { 'markersize':8,'marker':'h','label':'JEFF-3.2' },
					 { 'markersize':8,'marker':'d','label':'JEFF-3.2 + mod. Mugabghab' }
                     ]
		
df_measurements = pd.read_csv('Measurements.csv')
cleanupIsoNames(df_measurements)
df_measurements.sort_values(["Assembly","Sample","Isotope"])

#================================================
# TMI-1 NJ05YU 
#================================================
		
tmiTitle = "TMI-1 NJ05YU"
euFilename = "TMI1_Eu_CE"
gdFilename = "TMI1_Gd_CE"
df_TMI1 = pd.read_csv('TMI1_inventories.csv')
cleanupIsoNames(df_TMI1)

TMI_measured = df_measurements[df_measurements["Assembly"] == "TMI-1 NJ05YU"]
calculateCE(df_TMI1, TMI_measured, "mg / gIU")

TMI_burnups = {'AG536-C2D1': 53.5, 'AG536-C2D2':52.7, 'AG616-A1':45.9, 'AG616-B1':55.0, 'AG616-B2':52.4}
df_TMI1["Burnup"] = df_TMI1["Sample"].map(TMI_burnups)

isEu_TMI = df_TMI1.loc[(df_TMI1["Element"] == "Eu") & (df_TMI1["Isotope"] != "Eu-152") & (df_TMI1["Evaluation"] != "JEFF-3.2")]
isGd_TMI = df_TMI1.loc[(df_TMI1["Element"] == "Gd") & (df_TMI1["Isotope"] != "Gd-152") & (df_TMI1["Isotope"] != "Gd-157") & (df_TMI1["Evaluation"] != "JEFF-3.2")]

#pd.set_option('display.max_rows', 500)
#print(df_TMI1[["Sample","Evaluation","Isotope","Element","C/E"]])

plotIsosFP(isEu_TMI, tmiTitle, defaultSeriesInfo, euFilename, TMI_burnups, plotAspect=1.25,colWrap=3)
plotIsosFP(isGd_TMI, tmiTitle, defaultSeriesInfo, gdFilename, TMI_burnups, plotAspect=1.25,colWrap=3)
generateLatexTables(df_TMI1, TMI_measured, "TMI1_NJ05YU")

#================================================
# Calvert Cliffs
#================================================

ccTitle = "Calvert Cliffs"
euFilename = "CC_Eu_CE"
gdFilename = "CC_Gd_CE"
df_CC = pd.read_csv('CC_inventories.csv')
cleanupIsoNames(df_CC)

CC_measured = df_measurements[df_measurements["Assembly"] == "Calvert Cliffs"]
calculateCE(df_CC, CC_measured, "g / g Uinit")

CC_burnups = {"87-63": 44.34, "87-72":37.12, "87-81":27.35}
df_CC["Burnup"] = df_CC["Sample"].map(CC_burnups)

isEu_CC = df_CC.loc[(df_CC["Element"] == "Eu") & (df_CC["Isotope"] != "Eu-152") & (df_CC["Isotope"] != "Eu-151") & (df_CC["Evaluation"] != "JEFF-3.2")]
isGd_CC = df_CC.loc[(df_CC["Element"] == "Gd") & (df_CC["Isotope"] != "Gd-160") & (df_CC["Evaluation"] != "JEFF-3.2")]

plotIsosFP(isEu_CC, ccTitle, defaultSeriesInfo, euFilename, CC_burnups, plotAspect=1.25)
plotIsosFP(isGd_CC, ccTitle, defaultSeriesInfo, gdFilename, CC_burnups, plotAspect=1.25)
generateLatexTables(df_CC, CC_measured, "CalvertCliffs")
		
#================================================
# REBUS
#================================================

rebusTitle = "REBUS"
euFilename = "REBUS_Eu_CE"
df_REBUS = pd.read_csv('REBUS_inventories.csv')
cleanupIsoNames(df_REBUS)

REBUS_measured = df_measurements[df_measurements["Assembly"] == "REBUS"]
calculateCE(df_REBUS, REBUS_measured, "g / gUI")

REBUS_burnups = {"GKN-II":59.656}
df_REBUS["Burnup"] = df_REBUS["Sample"].map(REBUS_burnups)

isEu_REBUS = df_REBUS.loc[(df_REBUS["Element"] == "Eu") & (df_REBUS["Isotope"] != "Eu-152") & (df_REBUS["Isotope"] != "Eu-151") & (df_REBUS["Evaluation"] != "JEFF-3.2")]

plotIsosFP(isEu_REBUS, rebusTitle, defaultSeriesInfo, euFilename, REBUS_burnups, fontScale=1.25)
generateLatexTables(df_REBUS, REBUS_measured, "REBUS")
		
#================================================
# ARIANE
#================================================

ARIANETitle = "ARIANE"
euFilename = "ARIANE_Eu_CE"
df_ARIANE = pd.read_csv('ARIANE_inventories.csv')
cleanupIsoNames(df_ARIANE)

ARIANE_measured = df_measurements[df_measurements["Assembly"] == "ARIANE"]
calculateCE(df_ARIANE, ARIANE_measured, "g / gUI")

ARIANE_burnups = {"GU1":53.331}
df_ARIANE["Burnup"] = df_ARIANE["Sample"].map(ARIANE_burnups)

isEu_ARIANE = df_ARIANE.loc[(df_ARIANE["Element"] == "Eu") & (df_ARIANE["Isotope"] != "Eu-152") & (df_ARIANE["Isotope"] != "Eu-151") & (df_ARIANE["Evaluation"] != "JEFF-3.2")]

plotIsosFP(isEu_ARIANE, ARIANETitle, defaultSeriesInfo, euFilename, ARIANE_burnups, fontScale=1.25)
generateLatexTables(df_ARIANE, ARIANE_measured, "ARIANE")


#================================================
# CANDU 28-element / Pickering-A 19558C
#================================================

CANDUTitle = "Pickering-A 19558C (CANDU 28-element)"
euFilename = "CANDU_Eu_CE"
df_CANDU = pd.read_csv('CANDU_inventories.csv')
cleanupIsoNames(df_CANDU)

CANDU_measured = df_measurements[df_measurements["Assembly"] == "CANDU28"]
calculateCE(df_CANDU, CANDU_measured, "Calculated")

CANDU_burnups = {"Pickering-A 19558C":9.208}
df_CANDU["Burnup"] = df_CANDU["Sample"].map(CANDU_burnups)

isEu_CANDU = df_CANDU.loc[(df_CANDU["Element"] == "Eu") & (df_CANDU["Isotope"] != "Eu-152") & (df_CANDU["Isotope"] != "Eu-151") & (df_CANDU["Evaluation"] != "JEFF-3.2")]

plotIsosFP(isEu_CANDU, CANDUTitle, defaultSeriesInfo, euFilename, CANDU_burnups, fontScale=1.25)
generateLatexTables(df_CANDU, CANDU_measured, "CANDU")

#plotIsos(isEu, labels, defaultSeriesInfo, CANDUTitle, euFilename, CANDU_burnups) 

#================================================
# Burnup trend analysis
#================================================

pd.set_option('display.max_rows', 500)

df_Master = pd.concat([df_TMI1, df_CC, df_REBUS, df_ARIANE], ignore_index=True)
df_Master_Eu153 = df_Master.loc[(df_Master["Isotope"] == "Eu-153") & ((df_Master["Evaluation"] == "ENDF VII.0") | (df_Master["Evaluation"] == "ENDF VII.1"))]
df_Master_Eu154 = df_Master.loc[(df_Master["Isotope"] == "Eu-154") & (df_Master["Evaluation"] != "JEFF-3.2")]
df_Master_Eu155 = df_Master.loc[(df_Master["Isotope"] == "Eu-155") & (df_Master["Evaluation"] != "JEFF-3.2")]
df_Master_Gd154 = df_Master.loc[(df_Master["Isotope"] == "Gd-154") & (df_Master["Evaluation"] != "JEFF-3.2")]
df_Master_Gd155 = df_Master.loc[(df_Master["Isotope"] == "Gd-155") & (df_Master["Evaluation"] != "JEFF-3.2")]


labels_Eu153 = labels[0:2]
#["ENDF VII.0","ENDF VII.1","ENDF VII.0 + mod-Wright","ENDF VII.0 + mod-Mughabghab", "ENDF VII.1 + mod-Mughabghab", "JEFF-3.2 + mod-Mughabghab"]

plotIsoBuFP(df_Master_Eu153, labels_Eu153, defaultSeriesInfo, "Eu-153",plotAspect=2.5)
plotIsoBuFP(df_Master_Eu154, labels, defaultSeriesInfo, "Eu-154",plotAspect=2.5, plotFit=True)
plotIsoBuFP(df_Master_Eu155, labels, defaultSeriesInfo, "Eu-155",plotAspect=2.5)
plotIsoBuFP(df_Master_Gd154, labels, defaultSeriesInfo, "Gd-154",plotAspect=2.5)
plotIsoBuFP(df_Master_Gd155, labels, defaultSeriesInfo, "Gd-155",plotAspect=2.5)
#print(df_EuMaster[["Sample","Isotope","Evaluation","Burnup","C/E"]])