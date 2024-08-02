import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import seaborn as sns
from matplotlib.ticker import NullFormatter
from scipy.signal import savgol_filter
from scipy import stats
from sklearn.linear_model import LinearRegression
import datetime
import imageio as iio

def SetrcParams(rcParams):
     # Set fonts
     rcParams["font.family"] = "serif"
     # rcParams["text.usetex"] = True ## turn of if LaTeX not installed
     rcParams["font.size"] = 10

     return rcParams
plt.rcParams = SetrcParams(plt.rcParams)##Set "particular" style

handleList = ['Table 2 Raw data rep 1', 'Table 3 Raw data rep 2', 'Table 4 Raw data rep 3']
colorslist = ['tab:blue', 'tab:red', 'tab:orange']
c=0

for handle in handleList:
	color=colorslist[c]
	print(color)
	replabel = "Replicate " + handle[4]

	# read and prepare data file
	df = pd.read_excel('Data/Dataset_S1.xlsx', sheet_name=handle, skiprows=2).set_index('Well').transpose().reset_index()
	df = df.iloc[:,1:]
	df = df[['A1', 'A2', 'A3', 'A4', 'A5', 'A7', 'A8', 'A9', 'A10', 'A11',
	     'B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B8', 'B9', 'B10', 'B11', 
	     'C1', 'C2', 'C3', 'C4', 'C5',  'C7', 'C8', 'C9', 'C10', 'C11',
	    'D1', 'D2', 'D3', 'D4', 'D5',  'D7', 'D8', 'D9', 'D10', 'D11', 
	    'E1', 'E2', 'E3', 'E4', 'E5', 'E7', 'E8', 'E9', 'E10', 'E11', 
	    'F1', 'F2', 'F3', 'F4', 'F5',  'F7', 'F8', 'F9', 'F10', 'F11', 
	    'G1', 'G2', 'G3', 'G4', 'G5', 'G7', 'G8', 'G9', 'G10', 'G11']]
	#filter any negative values
	df = df.clip(lower=0)
	#log transformation
	df = np.log(df)

	# read mapping, and interval files
	metadf = pd.read_excel('Data/Dataset_S1.xlsx', sheet_name='Table 5 Nutrient conc. map', skiprows=2).set_index('Well')
	intervalsall = pd.read_excel('Data/Dataset_S1.xlsx', sheet_name='Table 1 Growth rates', skiprows=11).set_index("Well")
	rephandle = handle[-1]

	# Initiate figure
	fig, ax = plt.subplots(7,10,figsize=(10,7), squeeze=True)
	filename = 'data/' + handle

	# Iterate through wells/subplots
	n=0

	newdict = {}

	for axis in ax.flatten():
	    # gather nutrient data
	    wellname = df.iloc[:,n].name
	    glucose = metadf.loc[wellname]['Glucose (mM)']
	    ammonium = metadf.loc[wellname]['Ammonium (mM)']

	    # Plot the log transformed data as a blue line
	    x3 = df.index.values
	    y3 = df.iloc[:,n]
	    axis.plot(x3,y3, label=replabel, color=color)

	    # Fit the growth rates unless the well is excluded from fitting
	    if handle == 'Table 2 Raw data rep 1':
	    	excludewells = ['A1','A2','A3','A4','A5','A7','A8','A9','A10','A11']
	    else:
	    	excludewells = ['A7','A8','A9','A10','A11']
	    if not wellname in excludewells:

	        intervalstart =intervalsall.loc[wellname]['Fitting interval start (reading number) rep ' + str(rephandle)]
	        intervalend = intervalsall.loc[wellname]['Fitting interval end (reading number) rep ' + str(rephandle)]

	        x4 = df.index.values[intervalstart: intervalend]
	        y4 = df[wellname][intervalstart: intervalend]
	        slope, intercept, r, p, se = stats.linregress(x4, y4)
	        growthratehour = slope*6
	        # write to dictionary for writing to file later
	        newdict[wellname] = [glucose, ammonium, r, growthratehour, intervalstart, intervalend]

	        # plot the manual growth rate as a black line
	        gry = slope*x4 + intercept
	        axis.plot(x4,gry, c='black')

	        # annotate with the wellname, growth rate r, and blue bar for the fitting interval
	        axis.axvspan(intervalstart, intervalend, alpha=0.5, color=color)
	        axis.annotate(round(r,3), xy=(50,-1))

	    # Making the plots pretty
	    
	    # Grey background for the the 0 nutrient wells
	    if "A" in df.columns[n] or "1" in df.columns[n] and "10" not in df.columns[n] and "11" not in df.columns[n]:
	        axis.set_facecolor("lightgray")
	    
	    # Ticks
	    axis.xaxis.set_major_formatter(NullFormatter())
	    axis.yaxis.set_major_formatter(NullFormatter())
	    axis.tick_params(bottom = False)
	    axis.set_xticks([])
	    axis.set_ylim(-2,11)
	    if "A" in wellname:
	        axis.set_title(glucose, fontsize=10)
	    if "1" in wellname and "10" not in wellname and "11" not in wellname:
	        axis.set_ylabel(ammonium, fontsize=10, rotation='horizontal')
	        axis.yaxis.set_label_coords(-0.45, 0.5)
	    if "G" in wellname:
	        axis.tick_params(bottom=True)
	        axis.set_xticks([0,42,84,126])
	        axis.set_xticklabels(np.array([0,7,14,21]), fontsize=8)
	    else:
	        axis.tick_params(bottom=True)
	        axis.set_xticks([0,42,84, 126])
	        axis.set_xticklabels([], fontsize=8)
	    if "11" in wellname:
	        axis.set_yticks([1,5,10])
	        axis.set_yticklabels([1,5,10], fontsize=6)
	        axis.yaxis.set_label_position("right")
	        axis.yaxis.tick_right()
	    else:
	        axis.set_yticks([1,5,10])
	        axis.yaxis.tick_right()
	        axis.yaxis.set_label_position("right")
	        axis.set_yticklabels([], fontsize=6)
	        
	    n+=1

	# Make figure pretty
	fig.text(0.06, 0.5, 'Added [ammonium] (uM)', va='center', rotation='vertical', style='italic')
	fig.text(0.5, 0.92, 'Added [glucose] (uM)', ha='center', style='italic')    
	fig.text(0.5, 0.02, 'Time (hours)', ha='center')
	fig.text(0.98, 0.5, 'Log Luminescence value', va='center', rotation='vertical')
	plt.legend(bbox_to_anchor=[1.2,-0.3], fontsize=8)

	# Save things
	savehandle = handle[:-4] + 'finalintervals_plot.pdf'
	plt.legend().set_visible(False)
	plt.savefig(savehandle, dpi=400)
	c+=1
