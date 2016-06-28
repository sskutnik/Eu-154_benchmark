import matplotlib 
import matplotlib.pyplot as plt
import pandas as pd
#import numpy as np
import seaborn as sns
import re

sns.set(style="ticks",font_scale=1.5)

IMAGE_DEST='./images/'
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

	#print(df_calc[evalUnits])
	CE_series = []
	sigma_series = []

	# Match measurements to calculations by sample
	for sample in df_calc["Sample"].unique():
		df_measSample = df_measured[df_measured["Sample"] == sample].reset_index()
		# Loop over individual evaluations within samples for C/E values
		for eval in df_calc["Evaluation"].unique():
			isEval = df_calc[(df_calc["Evaluation"] == eval) & (df_calc["Sample"] == sample)].reset_index()
			tmpSeries = pd.Series(isEval[evalUnits] / df_measSample["Inventory"], name = "C/E")
			#print(isEval[evalUnits].count(), df_measSample.count())
			CE_series.append(tmpSeries)
			sigma_series.append(pd.Series(df_measSample["Uncertainty"], name = "Sigma"))
			
	# Add the C/E & measurement uncertainty data back to the measurements			
	joinedCE = pd.concat(CE_series, ignore_index=True)
	joinedSigma = pd.concat(sigma_series, ignore_index=True)
	df_calc["C/E"] = joinedCE
	df_calc["Sigma"] = joinedSigma 

	
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

	
def plotIsos(df_iso, labels, seriesInfo, pltTitle, filename, sigmaLabel="Sigma", yLim=DEFAULT_YLIM, legendLoc="best",legendAnchor=(1.1,0.5)):
	
	numIsos = len(df_iso["Isotope"].unique())
	isoIndex = range(0,len(df_iso["Isotope"].unique()))
	
	for idx in range(0,len(labels)):
		df_evalData = df_iso.loc[df_iso["Evaluation"] == labels[idx]]
		# Using .values attribute to handle cases where we exclude the first element (e.g., U-233), which plt.errorbar tries to access
		#plt.errorbar(df_evalData["Isotope"].index % numIsos, df_evalData["C/E"].values,df_evalData[sigmaLabel].values,ls='',**seriesInfo[idx])
		plt.errorbar(isoIndex, df_evalData["C/E"].values,df_evalData[sigmaLabel].values,ls='',**seriesInfo[idx])

	#print(df_iso.index,df_iso.index,df_iso["Isotope"])
	#plt.xticks(df_iso.index,df_iso["Isotope"])
	#plt.xlim(min(df_iso.index)-1,max(df_iso.index)+1)
	plt.xticks(isoIndex,df_iso["Isotope"].unique())
	plt.xlim(isoIndex[0]-0.5,isoIndex[-1]+2)
	plt.title(pltTitle)
	plt.ylabel("C/E")
	#plt.ylim(yLim)
	#plt.ylim(bottom=0.8)
	if(legendLoc == "best"):
		plt.legend(loc=legendLoc,title="Evaluation")
	else:
		plt.legend(loc=legendLoc,bbox_to_anchor=(legendAnchor),title="Evaluation")
	
	sns.despine(top=True)
	processPlot(plt, filename)

def plotIsosFP(df_iso, pltTitle, seriesInfo, filename):
	
	isoNumbers = []
	idx = 0
	for isotope in df_iso["Isotope"].unique():
		isoNumbers.append((isotope, idx))
		idx = idx + 1
	isoDict = dict(isoNumbers)
	
	isoMarkers = ['o', 'o', 's', '^', 'v']
	#print(isoDict)
	#df_iso["IsoNum"] = df_iso.loc("Isotope")
	#df_iso.replace({"IsoNum": isoDict}, inplace=True)
	df_iso["IsoNum"] = (df_iso.replace({"Isotope": isoDict}, inplace=False))["Isotope"]

	g = sns.factorplot(data=df_iso,kind='point',col='Sample',x='Isotope',y='C/E',hue='Evaluation',legend=False,scale=1.0,size=4,
	                   col_wrap=GLOB_COL_WRAP,linestyles='None',aspect=1.15, markers=isoMarkers)
	g.map(plt.errorbar,"IsoNum","C/E",'Sigma',data=df_iso,ls='',ms=0,zorder=0,elinewidth=2,color="gray")
	g.set_xlabels('Isotope').set_titles("{col_name}")
	#plt.legend(loc='lower left')
	#f = plt.figure()
	#f.subplots_adjust(right=0.8)
	#plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5),title="Evaluation")
	plt.legend(bbox_to_anchor=(1.05, 0.90), loc=2, borderaxespad=0., title="Evaluation")
	sns.despine(top=True)
	
	processPlot(plt, filename)
#	for ax, plotTitle in zip(g.axes.flat, titles):
#		ax.set_title("{:s} Sample {:s}".format(plotTitle, "{col_name}") )
	
	
#	plt.legend(loc=legendLoc)
def processPlot(plt, filename):
	plt.tight_layout()
	if SHOW_PLOTS:
		plt.show()
	else:
		plt.savefig(IMAGE_DEST + filename + '.' + OUTPUT_FORMAT)
		plt.close()
		
labels = ["ENDF VII.0","ENDF VII.0 + mod-Wright","ENDF VII.0 + mod-Mughabghab v1", "ENDF VII.1","ENDF VII.1 + mod-Mughabghab v1"]
# Default series meta-info
defaultSeriesInfo = [{ 'markersize':8,'marker':'o','label':'ENDF VII.0' },
                     { 'markersize':8,'marker':'o','label':'ENDF VII.1' },					 
                     { 'markersize':8,'marker':'^','label':'VII.0 + mod-Wright' },
                     { 'markersize':8,'marker':'v','label':'VII.0 + mod-Mughabghab v1' },
                     { 'markersize':8,'marker':'>','label':'VII.1 + mod-Mughabghab v1' },
                     ]
		
df_measurements = pd.read_csv('Measurements.csv')
cleanupIsoNames(df_measurements)
#df_measurements.sort_values(["Assembly","Sample","Isotope"])
#print(df_measurements[["Sample","Isotope","Element","Inventory"]])


#================================================
# TMI-1 NJ05YU 
#================================================
		
df_TMI1 = pd.read_csv('TMI1_inventories.csv')
cleanupIsoNames(df_TMI1)

TMI_measured = df_measurements[df_measurements["Assembly"] == "TMI-1 NJ05YU"]
calculateCE(df_TMI1, TMI_measured, "mg / gIU")

#pd.set_option('display.max_rows', 500)
#print(df_TMI1[["Sample","Isotope","Evaluation","Element","C/E"]])


isEu = df_TMI1.loc[(df_TMI1["Element"] == "Eu") & (df_TMI1["Isotope"] != "Eu-152")]
tmiTitle = "TMI-1 NJ05YU"
euFilename = "TMI1_Eu_CE"
plotIsosFP(isEu, tmiTitle, defaultSeriesInfo, euFilename)

#================================================
# Calvert Cliffs
#================================================

df_CC = pd.read_csv('CC inventories.csv')
cleanupIsoNames(df_CC)


CC_measured = df_measurements[df_measurements["Assembly"] == "Calvert Cliffs"]
calculateCE(df_CC, CC_measured, "g / g Uinit")

isEu = df_CC.loc[(df_CC["Element"] == "Eu") & (df_CC["Isotope"] != "Eu-152") & (df_CC["Isotope"] != "Eu-151")]
ccTitle = "Calvert Cliffs"
euFilename = "CC_Eu_CE"
plotIsosFP(isEu, ccTitle, defaultSeriesInfo, euFilename)
		
#================================================
# REBUS
#================================================

df_REBUS = pd.read_csv('REBUS_inventories.csv')
cleanupIsoNames(df_REBUS)

REBUS_measured = df_measurements[df_measurements["Assembly"] == "REBUS"]
calculateCE(df_REBUS, REBUS_measured, "g / gUI")

isEu = df_REBUS.loc[(df_REBUS["Element"] == "Eu") & (df_REBUS["Isotope"] != "Eu-152") & (df_REBUS["Isotope"] != "Eu-151")]
rebusTitle = "REBUS Sample GKN-II"
euFilename = "REBUS_Eu_CE"
#plotIsosFP(isEu, rebusTitle, defaultSeriesInfo, euFilename)
plotIsos(isEu, labels, defaultSeriesInfo, rebusTitle, euFilename) 

		
#================================================
# ARIANE
#================================================

df_ARIANE = pd.read_csv('ARIANE_inventories.csv')
cleanupIsoNames(df_ARIANE)


ARIANE_measured = df_measurements[df_measurements["Assembly"] == "ARIANE"]
calculateCE(df_ARIANE, ARIANE_measured, "g / gUI")
#(df_ARIANE["Isotope"] != "Eu-152") & (df_ARIANE["Isotope"] != "Eu-151")
isEu = df_ARIANE.loc[(df_ARIANE["Element"] == "Eu") & (df_ARIANE["Isotope"] != "Eu-152") & (df_ARIANE["Isotope"] != "Eu-151")]
ARIANETitle = "ARIANE Sample GU1"
euFilename = "ARIANE_Eu_CE"
#plotIsosFP(isEu, ARIANETitle, defaultSeriesInfo, euFilename)
plotIsos(isEu, labels, defaultSeriesInfo, ARIANETitle, euFilename) 

#================================================
# CANDU 28-element / Pickering-A 19558C
#================================================

df_CANDU = pd.read_csv('CANDU_inventories.csv')
cleanupIsoNames(df_CANDU)

CANDU_measured = df_measurements[df_measurements["Assembly"] == "CANDU28"]
calculateCE(df_CANDU, CANDU_measured, "Calculated")

#(df_CANDU["Isotope"] != "Eu-152") & (df_CANDU["Isotope"] != "Eu-151")
isEu = df_CANDU.loc[(df_CANDU["Element"] == "Eu") & (df_CANDU["Isotope"] != "Eu-152") & (df_CANDU["Isotope"] != "Eu-151")]
CANDUTitle = "Pickering-A 19558C (CANDU 28-element)"
euFilename = "CANDU_Eu_CE"
#plotIsosFP(isEu, CANDUTitle, defaultSeriesInfo, euFilename, legendLoc="right")
plotIsos(isEu, labels, defaultSeriesInfo, CANDUTitle, euFilename) 