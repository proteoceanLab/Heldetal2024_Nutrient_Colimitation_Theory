import numpy
from matplotlib import pyplot
import matplotlib
import colimitation_models
import colimitation_data_analysis
import colimitation_plots
import pandas


def main():

     # Number of replicate experiments
     num_replicates = 3
     rep_colors = ["tab:blue", "tab:red", "tab:orange"]

     # Read virtual supplementations for each data replicate 
     lim_data = [pandas.read_csv("./Data/Yield_data_lim_coeffs_rep" + str(r + 1) + ".csv", index_col=False) for r in range(num_replicates)]

     # L_threshold scan over a larger range to test number of virtual supplementations
     # above that threshold in experimental data
     L_threshold_range_large = numpy.linspace(0, 0.8, 17)
     num_abs_colim_data = numpy.array([[((lim_data[r]["L1"] > L_threshold) & (lim_data[r]["L2"] > L_threshold)).sum() for r in range(num_replicates)] for L_threshold in L_threshold_range_large])
 
     # L_threshold scan over a smaller range (based on the larger range) to test
     # fraction of virtual supplementations above that threshold in both 
     # experimental data and a simulated null model without colimitation
     L_threshold_range_small = numpy.linspace(0, 0.2, 11)
     fracs_abs_colim_data = numpy.array([[((lim_data[r]["L1"] > L_threshold) & (lim_data[r]["L2"] > L_threshold)).sum()/len(lim_data[r]["Meff"]) for r in range(num_replicates)] for L_threshold in L_threshold_range_small])

     # Read simulation data and calculate p-values
     num_sims = 10000
     null_model = "liebig_monod"
     pvalues_abs_colim = [[] for r in range(num_replicates)]
     fracs_abs_colim_sims_list = []

     # Read fractions of virtual supplementations for all thresholds
     dataframe = pandas.read_csv(f"./Data/Yield_sim_{null_model}_colim_fracs.csv", index_col=False)

     # Iterate over L_threshold
     for i in range(len(L_threshold_range_small)):
          L_threshold = L_threshold_range_small[i]
        
          # Get fractions for this threshold
          fracs_abs_colim_sims = list(dataframe[f"Fraction of virtual supplementations with L1 and L2 > {L_threshold}"])
          fracs_abs_colim_sims_list.append(fracs_abs_colim_sims)

          # Now calculate p-value for each experimental replicate against this
          # null model
          for r in range(num_replicates):
               pvalue = sum(int(f > fracs_abs_colim_data[i][r]) for f in fracs_abs_colim_sims)/len(fracs_abs_colim_sims)
               pvalues_abs_colim[r].append(pvalue)
     pvalues_abs_colim = numpy.array(pvalues_abs_colim)

################################################################################

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Initialize figure
     figure = pyplot.figure(figsize=(4*3, 3))
     figure.subplots_adjust(wspace=0.5)

     # Plot number of virtual supplementations with L1 and L2 both above a threshold
     axis = figure.add_subplot(1, 3, 1)
     axis.text(-0.25, 1.1, "A", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Threshold yield limitation coefficient $L^\mathrm{yield}$", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("Number of supplementations with\n$L^\mathrm{yield}_\mathrm{glu}$ and $L^\mathrm{yield}_\mathrm{amm}$ above threshold", fontsize=colimitation_plots.axis_label_size)
     axis.set_yscale("symlog", linthresh=1)
     axis.set_xlim([min(L_threshold_range_large), max(L_threshold_range_large)])
     axis.set_ylim([0, 1e3])
     for r in range(num_replicates):
          axis.plot(L_threshold_range_large, num_abs_colim_data[:, r], "-o", markersize=3, color=rep_colors[r], label="Rep " + str(r + 1))
          axis.text(0.98, 0.98 - 0.07*r, r"\textbf{Rep " + str(r + 1) + "}", color=rep_colors[r], transform=axis.transAxes, ha="right", va="top", fontsize=colimitation_plots.legend_label_size)

     # Plot fraction of virtual supplementations with L1 and L2 both above a threshold
     # For both experimental replicates and simulations
     axis = figure.add_subplot(1, 3, 2)
     axis.text(-0.25, 1.1, "B", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Threshold yield limitation coefficient $L^\mathrm{yield}$", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("Fraction of supplementations with\n$L^\mathrm{yield}_\mathrm{glu}$ and $L^\mathrm{yield}_\mathrm{amm}$ above threshold", fontsize=colimitation_plots.axis_label_size)
     axis.set_xlim([-0.01, 0.21])
     axis.set_ylim([0, 1])
     for r in range(num_replicates):
          axis.plot(L_threshold_range_small, fracs_abs_colim_data[:, r], "-o", markersize=3, color=rep_colors[r], label="Rep " + str(r + 1))
          axis.text(0.98, 0.98 - 0.07*r, r"\textbf{Rep " + str(r + 1) + "}", color=rep_colors[r], transform=axis.transAxes, ha="right", va="top", fontsize=colimitation_plots.legend_label_size)
     axis.text(0.98, 0.98 - 0.07*num_replicates, r"\textbf{Null model}" + "\n" + r"\textbf{no colim}", color="tab:purple", transform=axis.transAxes, ha="right", va="top", fontsize=colimitation_plots.legend_label_size)
     boxes = axis.boxplot(fracs_abs_colim_sims_list, positions=L_threshold_range_small, widths=0.012, sym="", manage_ticks=False, patch_artist=True, zorder=0, 
          medianprops={"color": "tab:purple", "linewidth": 0.5},
          boxprops={"facecolor": "tab:purple", "edgecolor": "tab:purple", "linewidth": 0.5},
          whiskerprops={"color": "tab:purple", "linewidth": 1.5},
          capprops={"color": "tab:purple", "linewidth": 1.5}, label="Null model\n(no colim)")
     axis.set_yscale("symlog", linthresh=1e-2)

     # Plot p-values for each experimental replicate against the null model at 
     # each threshold
     axis = figure.add_subplot(1, 3, 3)
     axis.text(-0.25, 1.1, "C", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Threshold yield limitation coefficient $L^\mathrm{yield}$", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("$p$-value compared to null model (no colim)", fontsize=colimitation_plots.axis_label_size)
     axis.set_yscale("symlog", linthresh=1e-3)
     axis.set_ylim([0, 1])
     axis.set_xlim([min(L_threshold_range_small), max(L_threshold_range_small)])
     for r in range(num_replicates):
          axis.plot(L_threshold_range_small, pvalues_abs_colim[r], "-o", markersize=3, color=rep_colors[r], label="Rep " + str(r + 1))  
          axis.text(0.98, 0.02 + 0.07*(num_replicates - r), r"\textbf{Rep " + str(r + 1) + "}", color=rep_colors[r], transform=axis.transAxes, ha="right", va="top", fontsize=colimitation_plots.legend_label_size)     

     figure.savefig("Figure_S30.pdf", bbox_inches="tight")

#################################################################################


if __name__ == '__main__':
     main()
