import colimitation_data_analysis
import colimitation_models
import colimitation_plots
import matplotlib
from matplotlib import pyplot
from matplotlib import colors
from matplotlib.gridspec import GridSpec
import numpy
import pandas


def main():

     # Number of replicate experiments
     num_replicates = 3
     rep_colors = ["tab:blue", "tab:red", "tab:orange"]

     # Number of R1 (glucose) and R2 (ammonium) values
     num_R1_values = 11
     num_R2_values = 7

     # Column labels for glucose and ammonium concentrations in data file
     R1_label = "Glucose (mM)"
     R2_label = "Ammonium (mM)"

     # Column labels for replicate yield measurements in data file
     yield_labels = ["Growth yield (OD 600 nm) rep " + str(r + 1) for r in range(num_replicates)]

     # Read R1, R2, and yield data into 2D meshes
     skiprows = 11
     sheet_name = "Table 1 Growth yields"
     R1_mesh_data, R2_mesh_data, yield_mesh_data = colimitation_data_analysis.Read2DScan("./Data/Dataset_S2.xlsx", R1_label, R2_label, yield_labels, num_R1_values, num_R2_values, sheet_name=sheet_name, skiprows=skiprows) 
     R1R2_list_data, yield_list_data = colimitation_data_analysis.Convert2DScanDataForFit(R1_mesh_data, R2_mesh_data, yield_mesh_data, reps=range(num_replicates))

     # Get parameters for a model fitted to the experimental data
     model_to_plot = "pat"
     dataframe = pandas.read_csv("./Data/Yield_data_fits_rep123.csv", index_col=False)
     dataframe_filtered = dataframe[dataframe["Model name"] == model_to_plot]
     row_index = dataframe_filtered.index[0]
     Nmax, s1, s2 = [float(x) for x in str(dataframe_filtered["Parameter values"][row_index]).split(";")]
     params = [Nmax, s1, s2]

     # Set up meshes for plots of fitted model
     R1_list_model = numpy.linspace(0, max(R1R2_list_data[0]), 50)
     R2_list_model = numpy.linspace(0, max(R1R2_list_data[1]), 50)
     R1_mesh_model, R2_mesh_model = numpy.meshgrid(R1_list_model, R2_list_model)
     yield_mesh_model = colimitation_models.CalcTraitForFit((R1_mesh_model, R2_mesh_model), params, model_to_plot)

     # Read virtual supplementations for each data replicate 
     lim_data = [pandas.read_csv("./Data/Yield_data_lim_coeffs_rep" + str(r + 1) + ".csv", index_col=False) for r in range(num_replicates)]

     # Calculate fraction of virtual supplementations with relative colimitation
     # (Meff > 1) in each replicate
     fracs_rel_colim_data = [(lim_data[r]["Meff"] > 1).sum()/len(lim_data[r]["Meff"]) for r in range(num_replicates)]

     # Read distributions of relative and absolute colimitation fractions from simulations
     model_lim_coeffs_no_colim = "liebig_monod"
     model_lim_coeffs_with_colim = "pat"
     dataframe = pandas.read_csv(f"./Data/Yield_sim_{model_lim_coeffs_no_colim}_colim_fracs.csv", index_col=False)
     fracs_rel_colim_sims_no_colim = list(dataframe["Fraction of virtual supplementations with L1 and L2 > 0.0"])
     dataframe = pandas.read_csv(f"./Data/Yield_sim_{model_lim_coeffs_with_colim}_colim_fracs.csv", index_col=False)
     fracs_rel_colim_sims_with_colim = list(dataframe["Fraction of virtual supplementations with L1 and L2 > 0.0"])   

     # Calculate p-values for relative and absolute colimitation for each data replicate
     pvalues_rel_colim = []
     print("p-values for fraction of virtual supplementations with colimitation")
     for r in range(num_replicates):
          pvalues_rel_colim.append(sum(int(f > fracs_rel_colim_data[r]) for f in fracs_rel_colim_sims_no_colim)/len(fracs_rel_colim_sims_no_colim))
          print("\t", "Rep " + str(r + 1) + ":", pvalues_rel_colim[-1])

################################################################################

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Initialize figure
     figure = pyplot.figure(figsize=(4*4, 3*1))
     figure.subplots_adjust(wspace=0.45, hspace=0.1)

     gs = GridSpec(2, 4, width_ratios=[1, 1, 1, 1], height_ratios=[0.05, 1], figure=figure)

################################################################################

     # Plot growth yield 2D scan data and fitted model
     axis = figure.add_subplot(gs[1, 0], projection="3d")
     axis.set_xlabel("Added glucose\nconcentration (mM)", fontsize=1.2*colimitation_plots.axis_label_size)
     axis.set_ylabel("Added ammonium\nconcentration (mM)", fontsize=1.2*colimitation_plots.axis_label_size)
     axis.set_zlabel("Growth yield\n(OD 600 nm)", fontsize=1.2*colimitation_plots.axis_label_size)
     axis.view_init(30, -135)
     axis.scatter(R1R2_list_data[0], R1R2_list_data[1], yield_list_data, color="tab:red", label="Data")
     axis.plot_surface(R1_mesh_model, R2_mesh_model, yield_mesh_model, label="Model fit")
     axis.legend(loc=(-0.3, 0.8), fontsize=1.2*colimitation_plots.legend_label_size)

################################################################################

     # Plot yield 1D scans
     axis = figure.add_subplot(gs[1, 1])
     axis.text(-0.25, 1.05, "B", transform=axis.transAxes, fontsize=1.2*colimitation_plots.panel_label_size)
     axis.set_xlabel("Added glucose concentration (mM)", fontsize=1.2*colimitation_plots.axis_label_size)
     axis.set_ylabel("Growth yield (OD 600 nm)", fontsize=1.2*colimitation_plots.axis_label_size)
     #axis.set_xticks([0, 1e-2, 1e-1, 1])
     axis.set_yticks(numpy.linspace(0, 0.6, 7))
     axis.set_ylim([-0.01, 0.6])

     norm = colors.SymLogNorm(linthresh=R2_mesh_data[:, 0][1], vmin=min(R2_mesh_data[:, 0]), vmax=max(R2_mesh_data[:, 0]))
     #norm = pyplot.Normalize(min(myvalues), max(myvalues))
     cmap = matplotlib.colormaps["viridis"]
     sm = pyplot.cm.ScalarMappable(cmap=cmap, norm=norm)
     sm.set_array([])
     ammonium_colors = sm.to_rgba(R2_mesh_data[:, 0])
     for i in range(len(R1_mesh_data)):
          #axis.errorbar(R1_mesh_data[i], yield_mesh_data_means[i], yerr=yield_mesh_data_sds[i], color=ammonium_colors[i], marker="o", markersize=5, linestyle="none")
          for r in range(num_replicates):
               xys = numpy.array([[R1_mesh_data[i][j], yield_mesh_data[i, j, r]] for j in range(len(R1_mesh_data[i])) if not isinstance(yield_mesh_data[i, j, r], str)])
               if len(xys) != 0:
                    axis.scatter(xys.T[0], xys.T[1], color=ammonium_colors[i], marker="o", s=5)
          R1s_model_fit = numpy.concatenate((numpy.linspace(R1_mesh_data[i][0], R1_mesh_data[i][1], 20), numpy.logspace(numpy.log10(R1_mesh_data[i][1]), numpy.log10(max(R1_mesh_data[i])), 50)))
          #R1s_model_fit = numpy.linspace(0, max(R1_mesh[i]), 50)
          yields_model_fit = colimitation_models.CalcTraitForFit((R1s_model_fit, R2_mesh_data[i][0]), params, model_to_plot)
          axis.plot(R1s_model_fit, yields_model_fit, "-", color=ammonium_colors[i])
     axis.set_xscale("symlog", linthresh=R1_mesh_data[0][2])

     # Add colorbar
     axis = figure.add_subplot(gs[0, 1])
     colorbar = figure.colorbar(sm, cax=axis, orientation="horizontal", location="top")
     colorbar.set_label("Added ammonium concentration (mM)", fontsize=1.2*colimitation_plots.axis_label_size)

     # Insert panel label for first panel using this axis, so that it aligns and prevents figure from cutting off
     # (These are problems with the 3D plot in panel A)
     axis.text(-1.75, 1.05, "A", transform=axis.transAxes, fontsize=1.2*colimitation_plots.panel_label_size)

################################################################################

     # Plot L1-L2 scatter plots
     axis = figure.add_subplot(gs[1, 2])
     axis.text(-0.32, 1.05, "C", transform=axis.transAxes, fontsize=1.2*colimitation_plots.panel_label_size)
     axis.set_xlabel("Growth yield limitation coefficient\nfor glucose $L_\mathrm{glu}^\mathrm{yield}$", fontsize=1.2*colimitation_plots.axis_label_size)
     axis.set_ylabel("Growth yield limitation\ncoefficient for ammonium $L_\mathrm{amm}^\mathrm{yield}$", fontsize=1.2*colimitation_plots.axis_label_size, labelpad=0)
     for r in range(num_replicates):
          axis.scatter(lim_data[r]["L1"], lim_data[r]["L2"], s=10, label="Rep " + str(r + 1), color=rep_colors[r])
          axis.text(0.02, 0.98 - 0.08*r, r"\textbf{Rep " + str(r + 1) + "}", color=rep_colors[r], transform=axis.transAxes, ha="left", va="top", fontsize=1.2*colimitation_plots.legend_label_size)
     axis.set_xlim([-1, 2.5])
     axis.set_ylim([-0.6, 1.2])
     axis.fill_between([0, 2.5], [0, 0], y2=1.2, color="0.9", zorder=-1)
     axis.text(1.25, 1, r"\textbf{Colimitation}", fontsize=1.2*colimitation_plots.axis_label_size, ha="center", va="top")

################################################################################

     # Plot Meff distributions
     axis = figure.add_subplot(gs[1, 3])
     axis.text(-0.28, 1.05, "D", transform=axis.transAxes, fontsize=1.2*colimitation_plots.panel_label_size)
     axis.set_xlabel("Fraction of virtual supplementations\nwith glucose-ammonium\ncolimitation ($L^\mathrm{yield}_\mathrm{glu},L^\mathrm{yield}_\mathrm{amm} > 0$)", fontsize=1.2*colimitation_plots.axis_label_size)
     axis.set_ylabel("Probability density", fontsize=1.2*colimitation_plots.axis_label_size)
     for r in range(num_replicates):
          axis.axvline(fracs_rel_colim_data[r], linestyle="--", color=rep_colors[r])
     axis.set_yscale("log")
     axis.set_xlim([0, 0.8])
     axis.set_ylim([5e-3, 1e1])
     axis.text(0.02, 0.98, r"\textbf{Sims}" + "\n" + r"\textbf{w/o colim}", color="tab:purple", transform=axis.transAxes, ha="left", va="top", fontsize=1.2*colimitation_plots.legend_label_size)
     axis.text(0.02, 0.83, r"\textbf{Sims}" + "\n" + r"\textbf{w/ colim}", color="tab:green", transform=axis.transAxes, ha="left", va="top", fontsize=1.2*colimitation_plots.legend_label_size)

     axis.hist(fracs_rel_colim_sims_no_colim, bins=numpy.arange(*axis.get_xlim(), 0.05), alpha=0.5, density=True, color="tab:purple", label="Sims\nw/o colim)")
     axis.hist(fracs_rel_colim_sims_with_colim, bins=numpy.arange(*axis.get_xlim(), 0.05), alpha=0.5, density=True, color="tab:green", label="Sims\nw/ colim)")
     axis.text(fracs_rel_colim_data[0] + 0.01, 1.1*axis.get_ylim()[1], r"\textbf{Rep 1}", color=rep_colors[0], rotation=90, ha="left", va="bottom", fontsize=1.2*colimitation_plots.axis_label_size)
     axis.text(fracs_rel_colim_data[1] - 0.01, 1.1*axis.get_ylim()[1], r"\textbf{Rep 2}", color=rep_colors[1], rotation=90, ha="right", va="bottom", fontsize=1.2*colimitation_plots.axis_label_size)
     axis.text(fracs_rel_colim_data[2] + 0.01, 1.1*axis.get_ylim()[1], r"\textbf{Rep 3}", color=rep_colors[2], rotation=90, ha="left", va="bottom", fontsize=1.2*colimitation_plots.axis_label_size)

################################################################################

     figure.savefig("Figure_S20.pdf", bbox_inches="tight")

################################################################################


if __name__ == '__main__':
     main()
