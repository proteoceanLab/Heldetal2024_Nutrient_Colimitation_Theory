import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import colimitation
import scipy
import numpy
from matplotlib.lines import Line2D
from matplotlib import gridspec
import matplotlib.font_manager as font_manager
import matplotlib.lines as mlines
from matplotlib import ticker as mticker
from matplotlib.ticker import StrMethodFormatter, NullFormatter
from scipy.stats import pearsonr, spearmanr, linregress
import colimitation_plots
from matplotlib.ticker import ScalarFormatter


plt.rcParams = colimitation_plots.SetrcParams(plt.rcParams)
'''plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 20})'''
plt.rcParams['lines.linewidth'] = 2


df = pd.read_excel("Dataset_S3.xlsx", skiprows = 11)
df = df.iloc[:-2,]

# fill nas and convert to uM for plotting

df['Nutrient 1 range (mM)'].fillna(0, inplace=True)
df['Nutrient 2 range (mM)'].fillna(0, inplace=True)

#convert to micromolar
df['Nutrient 1 K common unit (mM)'] = df['Nutrient 1 K common unit (mM)']*1000
df['Nutrient 2 K common unit (mM)'] = df['Nutrient 2 K common unit (mM)']*1000
df['Nutrient 1 environment In situ conc (mM)'] = df['Nutrient 1 environment In situ conc (mM)']*1000
df['Nutrient 2 environment In situ conc (mM)'] = df['Nutrient 2 environment In situ conc (mM)']*1000
df['Nutrient 1 range (mM)'] = df['Nutrient 1 range (mM)']*1000
df['Nutrient 2 range (mM)'] = df['Nutrient 2 range (mM)']*1000

plt.rcParams.update({'font.size': 12})
plt.rcParams['legend.fontsize'] = 9

fig = plt.figure(figsize=(8.2,4.5)) 
gs = gridspec.GridSpec(2,2) 
ax0 = plt.subplot(gs[0:2,0])
ax1 = plt.subplot(gs[0,1])
ax2 = plt.subplot(gs[1,1])

axes = [ax0, ax1, ax2]

# Stack on species and limiting nutrients
combinecolslist = [['Limiting Nutrient 1', 'Limiting Nutrient 2'],['Nutrient 1 K common unit (mM)','Nutrient 2 K common unit (mM)'],['Nutrient 1 environment In situ conc (mM)','Nutrient 2 environment In situ conc (mM)'],
['Nutrient 1 conc low (mM)','Nutrient 2 conc low (mM)'],['Nutrient 1 conc hi (mM)','Nutrient 2 conc hi (mM)'],['Habitat','Habitat'],['Nutrient 1 Li_liebig','Nutrient 2 Li_liebig'],
['Nutrient 1 Li_additive','Nutrient 2 Li_additive'],['Nutrient 1 Li_multiplicative','Nutrient 2 Li_multiplicative'],
['Nutrient 1 Li_pat','Nutrient 2 Li_pat'],['Nutrient 1 range (mM)','Nutrient 2 range (mM)']]
colslist = []
for colList in combinecolslist:
    cols = [df[col].squeeze() for col in colList]
    nutdf = pd.concat(cols, ignore_index=True)
    colslist.append(nutdf)
df2 = pd.concat(colslist, axis=1)
df2.columns = ['Limiting Nutrient', 'K common unit (uM)', 'in situ conc', 'lowconc', 'highconc', 'Habitat', 'Li_liebig', 'Li_additive', 'Li_multiplicative', 'Li_PAT', 'conc range']

# Subplot A

# setup colors
color_labels = df2['Limiting Nutrient'].unique()
marker_labels = df2['Habitat'].unique()
rgb_values = sns.color_palette().as_hex()

color_map = dict(zip(color_labels, rgb_values))
marker_values = ["o","^","s"]
marker_map = dict(zip(marker_labels, marker_values))
colors = df2['Limiting Nutrient'].map(color_map)
markers = df2['Habitat'].map(marker_map)

ogcolors = colors
ogmarkers = markers

#ax.scatter(df['K common unit (uM)'],df['conc low'], color=colors, yerr=df['conc low'])
#ax.errorbar(df['K common unit (uM)'],df['conc low'], yerr=df['conc low'], fmt="o", color=colors)
df2['conc range'].fillna(0, inplace=True)
df2.to_csv('test.csv')
for index, row in df2.iterrows():
    markers, caps, bars = axes[0].errorbar(x=row['in situ conc'],xerr=row['conc range']/2, y=row['K common unit (uM)'], markersize=5, fmt=marker_map[row['Habitat']], c=color_map[row['Limiting Nutrient']])
    [bar.set_alpha(0.4) for bar in bars]

axes[0].plot([0,1e2], [0,1e2], color='k', linestyle='--', linewidth=2)

axes[0].set_yscale('log')
axes[0].set_xscale('log')
axes[0].set_yticks([0.01,0.1,1,10,100])

#axes[0].set_ylim(1e-5,1e2)
axes[0].set_ylabel('Half-saturation concentration $K$ ($\mu$M)')
axes[0].set_xlabel('Resource concentration $R$ ($\mu$M)')

axes[0].annotate('more\nlimiting',xy=(0.5,40),
             horizontalalignment='left',fontsize=10)
axes[0].arrow(40,55,-20,0,head_width=5,head_length=5)
axes[0].arrow(80,55,50,0,head_width=5,head_length=50)
axes[0].annotate('less\nlimiting',xy=(220,40),
             horizontalalignment='left',fontsize=10)
fig.text(0.22,0.24,r'Spearman $\rho$ = 0.75', fontsize=11)
fig.text(0.22,0.2,r'p-value = $6\times10^{-8}$', fontsize=11)

s, p = spearmanr(df2['K common unit (uM)'],df2['in situ conc'])
#slope, intercept, r, p, se = linregress(df2['in situ conc'], df2['K common unit (uM)'])
print(s, p)


custom_lines = [Line2D([0], [0], color=rgb_values[0], lw=2),
                Line2D([0], [0], color=rgb_values[1], lw=2),
                Line2D([0], [0], color=rgb_values[2], lw=2),
                Line2D([0], [0], color=rgb_values[3], lw=2),
                Line2D([0], [0], color=rgb_values[4], lw=2),
                Line2D([0], [0], color=rgb_values[5], lw=2),
                Line2D([0], [0], color=rgb_values[6], lw=2),
                Line2D([0], [0], color=rgb_values[7], lw=2),
                mlines.Line2D([], [], color='black', marker='o', linestyle='None',
                          markersize=5, label='Marine'),
                mlines.Line2D([], [], color='black', marker='^', linestyle='None',
                          markersize=5, label='Freshwater'),
                mlines.Line2D([], [], color='black', marker='s', linestyle='None',
                          markersize=5, label='Enteric')]
color_labels2 = np.append(color_labels, ['Marine', 'Freshwater', 'Enteric'])
#axes[1].legend(custom_lines, color_labels2, prop={'size': 9}, bbox_to_anchor=(1.12, 0.5))
marine = mlines.Line2D([], [], color='black', marker='o', linestyle='None',
                          markersize=5, label='Marine')
df2['Klog'] = df2['K common unit (uM)'].apply(lambda x: np.log(x))
df2['conclog'] = df2['in situ conc'].apply(lambda x: np.log(x))


x = []
y = []
species = []
nutrient = []
indexes = []
habitats = []
#Meff = []
n=0

for index, row in df2.iterrows():
    x += ['Liebig', 'PAT', 'Add.', 'Mult.']
    y += [row['Li_liebig'], row['Li_PAT'], row['Li_additive'], row['Li_multiplicative']]
    #species += [row['Species 2'], row['Species 2'], row['Species 2'], row['Species 2']]
    nutrient += [row['Limiting Nutrient'], row['Limiting Nutrient'], row['Limiting Nutrient'], row['Limiting Nutrient']]
    habitats += [row['Habitat'], row['Habitat'], row['Habitat'], row['Habitat']]
    #Meff += [row['Meff_liebigsmooth'], row['Meff_pat'], row['Meff_sequential'], row['Meff_multiplicative']]

jitterdf = pd.DataFrame(x, y)
#jitterdf['species'] = species
jitterdf['nutrient'] = nutrient
jitterdf['habitat'] = habitats
jitterdf.reset_index(inplace=True)
jitterdf.to_csv('jitter.csv')
jitterdf.columns = ['Li', 'model', 'nutrient', 'habitat']
jitterdf['marker'] = jitterdf['habitat'].map(marker_map)
#df_unit = jitterdf.groupby('model').median().reset_index()
#jitterdf = jitterdf.sort_values(by='habitat', ascending=True)
jitterdfgrouped = jitterdf.groupby('nutrient')
#colors = ['tab:red', 'tab:orange', 'tab:cyan']
n = 0
#blackbarsns.scatterplot(x="model", y="Li", data=df_unit, marker='_',s=600,color='k', ax=axes[1])
for name, group in jitterdfgrouped:
    marker = marker_map[group['habitat'].iloc[0]]
    color = color_map[name]
    sns.stripplot(data =group, x="model", y="Li", jitter=0.2, size=4, marker = marker, color=color, ax=axes[1])
    n += 1

axes[1].set_yscale('symlog', linthresh=0.00001)
axes[1].set_xlabel("")
axes[1].set_ylabel('Rate limitation\ncoefficient $L^\mathrm{rate}_i$')
plt.savefig("Figure3.pdf", dpi=400, bbox_inches='tight')

axes[1].set_ylim([0,1.1])
axes[1].set_yticks([0,0.00001,0.0001,0.001,0.01,0.1,1])
axes[1].legend(custom_lines, color_labels2, prop={'size': 9}, bbox_to_anchor=(1.1,1))

#axes[1].legend().set_visible(False)
custom_lines2 = [mlines.Line2D([], [], color='black', marker='o', linestyle='None',
                          markersize=5, label='Marine'),
                mlines.Line2D([], [], color='black', marker='^', linestyle='None',
                          markersize=5, label='Freshwater'),
                mlines.Line2D([], [], color='black', marker='s', linestyle='None',
                          markersize=5, label='Enteric')]
habitat_labels = ['Marine', 'Freshwater', 'Enteric']

Meffdf = df[['Meff_liebig', 'Meff_PAT','Meff_additive', 'Meff_multiplicative', 'Habitat']]
Meffdf.columns = ['Liebig', 'PAT', 'Add.', 'Mult.', 'habitat']
Meffdf = Meffdf.set_index("habitat")
Meffdf = Meffdf.stack().reset_index()

Meffdf.columns = ['habitat', 'model', 'Meff']
Meffdf['marker'] = Meffdf['habitat'].map(marker_map)

Meffgrouped = Meffdf.groupby('habitat')
n = 0
for name, group in Meffgrouped:
    marker = group['marker'].iloc[0]
    color = colors[n]
    sns.stripplot(data =group, x="model", y="Meff", jitter=0.2, size=4, marker = marker, color='black', ax=axes[2], alpha=0.5)
    n += 1
axes[2].set_xlabel("")
axes[2].set_ylabel("Number rate-limiting\nresources $M^\mathrm{rate}_\mathrm{eff}$")

fig.text(0.01,0.95,'A', fontsize=16)
fig.text(0.44,0.95,'B', fontsize=16)
fig.text(0.44,0.53,'C', fontsize=16)

plt.tight_layout(pad=2.0)
plt.savefig("Figure4.pdf", dpi=400, bbox_inches='tight')