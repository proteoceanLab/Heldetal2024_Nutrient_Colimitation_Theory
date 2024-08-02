import pandas as pd
import matplotlib.pyplot as plt
import colimitation_plots

plt.rcParams = colimitation_plots.SetrcParams(plt.rcParams)
'''plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 20})'''
plt.rcParams['lines.linewidth'] = 2

df = pd.read_excel("Data/Dataset_S3.xlsx", skiprows=4)
df = df.iloc[:-3,]
df2 = pd.DataFrame()
df2['Liebig_Monod'] = pd.concat([df['Nutrient 1 Li_liebig'], df['Nutrient 2 Li_liebig']])
df2['Multiplicative_Monod'] = pd.concat([df['Nutrient 1 Li_multiplicative'], df['Nutrient 2 Li_multiplicative']])
df2['Additive'] = pd.concat([df['Nutrient 1 Li_additive'], df['Nutrient 2 Li_additive']])
df2['PAT'] = pd.concat([df['Nutrient 1 Li_pat'], df['Nutrient 2 Li_pat']])


fig, ax = plt.subplots(4,2, figsize=(6,6))
ax[0,0].hist(df['Meff_liebig'], bins=5, color='#0504aa', align='mid', rwidth=0.85)
ax[1,0].hist(df['Meff_PAT'], bins=5, color='#0504aa', align='left',rwidth=0.95)
ax[2,0].hist(df['Meff_additive'], bins=5, color='#0504aa', align='left',rwidth=0.9)
ax[3,0].hist(df['Meff_multiplicative'], bins=5, color='#0504aa', align='left', rwidth=0.85)
ax[3,0].set_xlabel('Number of rate limiting resources $M^\mathrm{rate}_\mathrm{eff}$')

for axis in [ax[0,0],ax[1,0],ax[2,0],ax[3,0]]:
    axis.set_xlim([0.9,2.2])
    axis.set_ylim([0,22])
    axis.set_xticks([1,1.25,1.5,1.75,2])
    axis.set_yticks([0,5,10,15,20])

ax[0,1].hist(df2['Liebig_Monod'], bins=5, color='#0504aa', align='left')
ax[1,1].hist(df2['PAT'], bins=5, color='#0504aa', align='left',rwidth=0.8)
ax[2,1].hist(df2['Additive'], bins=5, color='#0504aa', align='left',rwidth=0.8)
ax[3,1].hist(df2['Multiplicative_Monod'], bins=5, color='#0504aa', align='left',rwidth=0.8)
ax[3,1].set_xlabel('Rate limitation coefficient $L^\mathrm{rate}i$')

rows = ['Liebig\nMonod', 'PAT', 'Additive', 'Mult.\nMonod']
for axis, row in zip(ax[:,1], rows):
    axis.set_ylabel(row, rotation=0)
    axis.yaxis.set_label_coords(1.2,0.4)

for axis in [ax[0,1],ax[1,1],ax[2,1],ax[3,1]]:
    axis.set_xlim([-0.1,1.1])
    axis.set_ylim([0,40])
    axis.set_xticks([0,0.25,0.5,0.75,1])
    axis.set_yticks([0,10,20,30,40])
fig.supylabel('No. of organism-resource combinations')
plt.tight_layout()

plt.savefig("FigureSX_Histograms.pdf")
