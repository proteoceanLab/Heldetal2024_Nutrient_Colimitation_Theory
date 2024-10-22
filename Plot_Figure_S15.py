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

     # Number of R1 (glucose) and R2 (ammonium) values
     num_R1_values = 10
     num_R2_values = 7

     # Column labels for glucose and ammonium concentrations in data file
     R1_label = "Glucose (mM)"
     R2_label = "Ammonium (mM)"

     # Read virtual supplementations for two simulated models without noise
     lim_sim_liebig_monod = pandas.read_csv("./Data/Rate_sim_liebig_monod_rmin_noiseless_lim_coeffs.csv", index_col=False)
     lim_sim_gm = pandas.read_csv("./Data/Rate_sim_pat_rmin_noiseless_lim_coeffs.csv", index_col=False)

     # Read virtual supplementations for each data replicate 
     lim_data = [pandas.read_csv("./Data/Rate_data_lim_coeffs_rep" + str(r + 1) + ".csv", index_col=False) for r in range(num_replicates)]

     # Read noisy simulated data and calculate limitation coefficients for sim 1
     num_sims = 10000
     example_sim = 0
     trait_labels = [f"Growth rate (per hour) sim {s + 1}" for s in range(num_sims)]
     R1_mesh_sim, R2_mesh_sim, trait_mesh_sim  = colimitation_data_analysis.Read2DScan(f"./Data/Rate_sim_liebig_monod_rmin.csv", R1_label, R2_label, trait_labels, num_R1_values, num_R2_values) 
     lim_sim = colimitation_data_analysis.CalculateVirtualSupplementations(R1_mesh_sim, R2_mesh_sim, trait_mesh_sim[:, :, example_sim])

################################################################################

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Initialize figure
     figure = pyplot.figure(figsize=(4*2, 3))
     figure.subplots_adjust(wspace=0.5, hspace=0.5)

     # Plot Meff distributions for simulated models without noise
     axis = figure.add_subplot(1, 2, 1)
     axis.text(-0.2, 1.05, "A", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Effective number of rate-limiting resources\n$M_\mathrm{eff}^\mathrm{rate}$ (glucose and ammonium)", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("Density of supplementations", fontsize=colimitation_plots.axis_label_size)
     axis.set_xlim([1, 2])
     axis.set_ylim([0, 5])
     axis.hist(lim_sim_gm["Meff"], bins=numpy.linspace(1, 2, 6), histtype="stepfilled", color="tab:green", density=True, alpha=0.5)
     axis.hist(lim_sim_liebig_monod["Meff"], bins=numpy.linspace(1, 2, 6), histtype="stepfilled", color="tab:purple", density=True, alpha=0.5)
     axis.text(0.98, 0.98, r"\textbf{Sim no noise with colim}" + "\n" + r"\textbf{(PAT model)}", color="tab:green", transform=axis.transAxes, ha="right", va="top", fontsize=colimitation_plots.legend_label_size)
     axis.text(0.98, 0.81, r"\textbf{Sim no noise no colim}" + "\n" + r"\textbf{(Liebig Monod model)}", color="tab:purple", transform=axis.transAxes, ha="right", va="top", fontsize=colimitation_plots.legend_label_size)

     # Plot Meff distributions for real data and one simulation with noise
     axis = figure.add_subplot(1, 2, 2)
     axis.text(-0.2, 1.05, "B", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Effective number of rate-limiting resources\n$M_\mathrm{eff}^\mathrm{rate}$ (glucose and ammonium)", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("Density of supplementations", fontsize=colimitation_plots.axis_label_size)
     axis.set_xlim([-2, 2])
     axis.set_ylim([0, 2])
     for r in range(num_replicates):
          axis.hist(lim_data[r]["Meff"], bins=numpy.linspace(-2, 2, 21), histtype="stepfilled", color=rep_colors[r], density=True)
          axis.text(0.02, 0.98 - 0.07*r, r"\textbf{Rep " + str(r + 1) + "}", color=rep_colors[r], transform=axis.transAxes, ha="left", va="top", fontsize=colimitation_plots.legend_label_size)
     axis.hist(lim_sim["Meff"], bins=numpy.linspace(-2, 2, 21), histtype="stepfilled", color="gray", density=True)
     axis.text(0.02, 0.98 - 0.07*num_replicates, r"\textbf{Sim w/ noise, no colim}" + "\n" + r"\textbf{(Liebig Monod model)}", color="gray", transform=axis.transAxes, ha="left", va="top", fontsize=colimitation_plots.legend_label_size)

     figure.savefig("Figure_S15.pdf", bbox_inches="tight")

#################################################################################


if __name__ == '__main__':
     main()