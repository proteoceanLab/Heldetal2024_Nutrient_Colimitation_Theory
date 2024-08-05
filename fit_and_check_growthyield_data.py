import pandas as pd
import numpy as np
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

mpl.rcParams.update(mpl.rcParamsDefault)

def SetrcParams(rcParams):
     # Set fonts
     rcParams["font.family"] = "serif"
     # rcParams["text.usetex"] = True ## turn of if LaTeX not installed
     rcParams["font.size"] = 10

     return rcParams
plt.rcParams = SetrcParams(plt.rcParams)##Set "particular" style

##Read in Data##
datafile="Data/Dataset_S2.xlsx"
rep1 = pd.read_excel(datafile, sheet_name="Table 2 Raw data rep 1", header=2,index_col=0).T
rep2 = pd.read_excel(datafile, sheet_name="Table 2 Raw data rep 1", header=2,index_col=0).T
rep3 = pd.read_excel(datafile, sheet_name="Table 2 Raw data rep 1", header=2,index_col=0).T
all_reps = [rep1, rep2, rep3]

conc_df= pd.read_excel("Data/Dataset_S2.xlsx", sheet_name='Table 5 Nutrient conc. map', skiprows=2).set_index('Well')
conc_df = conc_df.round(3)
##---------------------------------------------------------------------------##
##Calculate and compile data##
for df in all_reps:
    dt_index = pd.to_datetime(pd.Series(df.index),format="%H:%M:%S",errors="coerce")
    dt_index.fillna(pd.to_datetime("1900-01-01 23:59:59"),inplace=True)
    df.index = dt_index#.dt.time

regex = "\\D1$|H" #This will match any column that is called "x1", where x is any letter, or has the letter H
norm_reps = []
rep_ccs,rep_std,rep_devmean = [],[],[]
bg_factors = []
for i,rep in enumerate(all_reps):
    rep_bg = rep.filter(regex=regex)
    avg_bg = np.mean(rep_bg)
    norm_rep = rep - avg_bg

    carrying_cap = norm_rep.loc[dt.time(12,0,0):dt.time(16,0,0)].mean() #Take the average value for 12-16 hours to calculate carrying capactity
    carrying_capstdev = norm_rep.loc[dt.time(12,0,0):dt.time(16,0,0)].std()
    devmean = carrying_capstdev/carrying_cap

    carrying_cap[carrying_cap < 0] = 0

    calc_df = pd.concat([carrying_cap,carrying_capstdev,devmean],axis=1)
    calc_df.rename(columns={
            0:f"mean_biomass_yield_Rep-{i+1}",
            1:f"stdev_biomass_yield_Rep-{i+1}",
            2:f"stdev:mean_biomass_yield_Rep-{i+1}"
            },inplace=True
        )
    norm_reps.append(norm_rep)
    bg_factors.append(avg_bg)

    conc_df = pd.concat([conc_df,calc_df],axis=1)




conc_df.to_csv("processed_data.csv")##Save the calculated data
##---------------------------------------------------------------------------##

##Initialize plot##
fig, axs = plt.subplots(8,12,figsize=(12,8),squeeze=False) ##96 subplots
##---------------------------------------------------------------------------##
##Add labels and Concentration data to axes##
fig.text(0.09, 0.5, s="Ammonium added (mM)",va="center",rotation="vertical")
ammo_concs = sorted(conc_df["Ammonium (mM)"].unique().tolist())+[0.0]
for i,ammo in enumerate(ammo_concs):
    loc = np.linspace(.15,.85,8)
    fig.text(0.105, loc[i],s=ammo,va="center",rotation="vertical",fontsize=10)

fig.text(0.5, 0.92, s="Glucose added (mM)",ha="center")
gluc_concs = [0.0]+sorted(conc_df["Glucose (mM)"].unique().tolist(),reverse=True)
for i,gluc in enumerate(gluc_concs):
    loc = np.linspace(.15,.87,12)
    fig.text(loc[i],0.89,s=gluc,ha="center",fontsize=10)

fig.text(0.5, 0.05, s="Time (hours)",ha="center")
fig.text(0.92, 0.5, s="OD 600",va="center",rotation="vertical")
##---------------------------------------------------------------------------##
##Set look and colors
kwargs = {"fontsize":5,"rotation":90,"zorder":0,"va":"bottom"}
linestyle = {"c":"gray", "ls":"--","lw":1}
#colors = ("#1f77b4","#ff7f0e","#2ca02c")
colors = ("tab:blue","tab:red","tab:orange")
##---------------------------------------------------------------------------##
for type in ("linear", 'log'): ##Plot the data on linear and log axes
    for i,axis in enumerate(axs.flatten()):

        axis.set_yscale(type)
        if type  == "log":
            lims = [0.01,1]
            labels = [0.01,0.1,1]
        else:
            lims = [0,0.8]
            labels = [0,0.4,0.8]

        well = norm_reps[0].columns[i]
        ##Plot the data##
        for j in range(3):
            axis.plot(norm_reps[j][well],label=f"rep{j+1}", color=colors[j])
            axis.axhline(y = conc_df[f"mean_biomass_yield_Rep-{j+1}"][well],**linestyle)

        if (i < 12)|(i > 84)|(i % 12==0)|((i + 1) % 12==0):
            axis.set_facecolor("lightgray") ## make blank wells gray

        ## Format Y Axes on all subplots ##
        axis.set_ylim(lims)
        axis.yaxis.tick_right()
        axis.set_yticks(labels,labels=labels,**kwargs)

        axis.tick_params(axis="y", which="both",direction="in",length=4)
        axis.tick_params(axis="x", which="both",direction="in",length=2)

        axis.spines["left"].set_visible(False)
        axis.spines["top"].set_visible(False)

        ## Format X Axis on bottom subplots ##
        if i < 84:
            axis.set_xticklabels([])
        else:
            axis.xaxis.set_major_formatter(mdates.DateFormatter("%H"))
            xlabels = axis.xaxis.get_ticklabels()
            xlabels[-1].set_text("24")
            [xl.set_text("") for xl in xlabels[1::2]]
            axis.set_xticklabels(xlabels, fontsize=6,rotation=45)##Generates a UserWarning, can be ignored
    

    ##Create the linear and log versions of the graphs##
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    plt.savefig(f"FigureS13_{type}.pdf")
##---------------------------------------------------------------------------##
