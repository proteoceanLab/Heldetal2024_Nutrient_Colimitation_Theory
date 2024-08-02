import numpy
import matplotlib
from matplotlib import pyplot
from matplotlib import colors
import pandas
import colimitation_models
import colimitation_plots


def main():

     # Read data table from Moore et al. paper
     data = pandas.read_csv("./Data/Moore_etal_2013_NatureGeosci_Supplementary_Table_1.csv")

     # List of pairs: elements and relative biomass yields
     pair_list = []

     # Iterate over rows in data table
     for index, row in data.iterrows():

          # Skip elements with missing data for these entries (H and O)
          if (not numpy.isnan(row["Phytoplankton quota (mol:mol C)"])) and (not numpy.isnan(row["Mean ocean concentration (umol/kg)"])):

               # Get relative yield (biomass per unit element)
               Y = 1/row["Phytoplankton quota (mol:mol C)"]

               # Get element concentration and convert to mol from umol
               R = row["Mean ocean concentration (umol/kg)"]/1e6

               # Add data to list
               pair_list.append((row["Element"], R*Y))

     # Sort elements and ratios in order of increasing biomass
     pair_list = sorted(pair_list, key=lambda x: x[1])

     # Extract sorted lists of elements and yields separately
     elements = [pair[0] for pair in pair_list]
     yields = numpy.array([pair[1] for pair in pair_list]) #/pair_list[0][1]

     # Calculate Meff as a function of q
     recip_q_range = numpy.linspace(1e-1, 100, 100)
     Meff_range = [colimitation_models.CalcMeffGM(yields, numpy.ones(len(yields)), 1/recip_q) for recip_q in recip_q_range]

     # Calculate limitation coefficients for three values of q
     q1_index, q2_index, q3_index = 0, 2, len(recip_q_range) - 1
     Lis_q1 = [colimitation_models.CalcLimCoeffGM(i, yields, numpy.ones(len(yields)), 1/recip_q_range[q1_index]) for i in range(len(elements))]
     Lis_q2 = [colimitation_models.CalcLimCoeffGM(i, yields, numpy.ones(len(yields)), 1/recip_q_range[q2_index]) for i in range(len(elements))]
     Lis_q3 = [colimitation_models.CalcLimCoeffGM(i, yields, numpy.ones(len(yields)), 1/recip_q_range[q3_index]) for i in range(len(elements))] 

     # Set up log colormap for reciprocal q values
     norm = colors.LogNorm(vmin=min(recip_q_range), vmax=max(recip_q_range))
     cmap = matplotlib.colormaps["summer"]
     sm = pyplot.cm.ScalarMappable(cmap=cmap, norm=norm)
     sm.set_array([])
     qcolors = sm.to_rgba(recip_q_range)

################################################################################

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Initialize figure
     figure = pyplot.figure(figsize=(6, 3*2))
     figure.subplots_adjust(hspace=0.4)

     # Plot yields
     axis = figure.add_subplot(2, 1, 1)
     axis.text(-0.15, 1.1, "A", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.bar(range(len(yields)), yields, color="tab:purple")
     axis.set_xticks(range(len(yields)))
     axis.set_xticklabels(elements)
     axis.set_ylabel("Maximum potential yield\n(mol C per kg seawater)", fontsize=1.2*colimitation_plots.axis_label_size)
     axis.set_xlim([-0.5, len(yields) - 0.5])
     axis.set_ylim([1e-6, 1e2])
     axis.set_yscale("log")

     # Plot limitation coefficients
     axis = figure.add_subplot(2, 1, 2)
     axis.text(-0.15, 1.05, "B", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.bar(numpy.linspace(0, len(elements), len(elements)), Lis_q1, width=0.25, label="Weak\ninteraction", color=qcolors[q1_index])
     axis.bar(numpy.linspace(0, len(elements), len(elements)) + 0.3, Lis_q2, width=0.25, label="Intermediate\ninteraction", color=qcolors[q2_index])
     axis.bar(numpy.linspace(0, len(elements), len(elements)) + 0.6, Lis_q3, width=0.25, label="Strong\ninteraction", color=qcolors[q3_index])
     axis.set_xticks(numpy.linspace(0, len(elements), len(elements)) + 0.3)
     axis.set_xticklabels(elements)
     axis.set_ylabel("Yield limitation\ncoefficients $L^\mathrm{yield}_i$", fontsize=1.2*colimitation_plots.axis_label_size)
     axis.set_xlim([-0.4, len(yields) + 0.9])
     axis.legend(loc=(0.05, 0.4))

     # Plot Meff vs. q
     axis_inset = axis.inset_axes([0.5, 0.4, 0.48, 0.55])
     for i in range(len(recip_q_range) - 1):
          axis_inset.plot(recip_q_range[i:i + 2], Meff_range[i:i + 2], "-", color=qcolors[i], linewidth=5)
     axis_inset.set_xlabel("Resource interaction strength", fontsize=colimitation_plots.axis_label_size)
     axis_inset.set_ylabel("No. yield-limiting\nresources $M^\mathrm{yield}_\mathrm{eff}$", fontsize=colimitation_plots.axis_label_size)
     axis_inset.set_ylim([1, len(elements)])
     axis_inset.set_yticks(numpy.arange(1, len(elements) + 1, 2))
     axis_inset.set_yticklabels([str(y) for y in numpy.arange(1, len(elements) + 1, 2)])
     axis_inset.set_yticks(numpy.arange(1, len(elements) + 1, 2), minor=True)

     figure.savefig("Figure_5.pdf", bbox_inches="tight")

#################################################################################


if __name__ == '__main__':
     main()
