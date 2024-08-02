import colimitation_data_analysis
import colimitation_models
import colimitation_plots
import matplotlib
from matplotlib import pyplot
from matplotlib import colors
import numpy
import pandas


def main():

     # Number of replicate experiments
     num_replicates = 3

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

     # Calculate yield means and SDs over replicates
     yield_mesh_data_means = numpy.zeros(R1_mesh_data.shape)
     yield_mesh_data_sds = numpy.zeros(R1_mesh_data.shape)
     for i in range(len(yield_mesh_data)):
          for j in range(len(yield_mesh_data[i])):
               # Get all yield data for this well, skip any entries that are strings ("not fit")
               yields = [r for r in yield_mesh_data[i][j] if not isinstance(r, str)]
               if len(yields) != 0:
                    yield_mesh_data_means[i][j] = numpy.mean(yields)
                    yield_mesh_data_sds[i][j] = numpy.std(yields)
               else:
                    yield_mesh_data_means[i][j] = numpy.nan
                    yield_mesh_data_sds[i][j] = numpy.nan

################################################################################

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Initialize figure
     figure = pyplot.figure(figsize=(4*2, 3))
     figure.subplots_adjust(wspace=0.6)

     # Plot yield 1D scans over glucose
     axis = figure.add_subplot(1, 2, 1)
     axis.text(-0.3, 1.05, "A", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Reciprocal added glucose\nconcentration (1/mM)", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("Reciprocal growth yield (1/OD 600 nm)", fontsize=colimitation_plots.axis_label_size)
     axis.set_xlim([-0.4, 14])
     axis.set_xticks(numpy.linspace(0, 14, 8))
     axis.set_ylim([0, 80])

     norm = colors.SymLogNorm(linthresh=R2_mesh_data[:, 0][1], vmin=min(R2_mesh_data[:, 0]), vmax=max(R2_mesh_data[:, 0]))
     cmap = matplotlib.colormaps["viridis"]
     sm = pyplot.cm.ScalarMappable(cmap=cmap, norm=norm)
     sm.set_array([])
     ammonium_colors = sm.to_rgba(R2_mesh_data[:, 0])
     colorbar = figure.colorbar(sm, ax=axis)
     colorbar.set_label("Added ammonium\nconcentration (mM)", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=25)
     for i in range(len(yield_mesh_data_means)):
          #axis.errorbar(R1_mesh_data[i], yield_mesh_data_means[i], yerr=yield_mesh_data_sds[i], color=ammonium_colors[i], marker="o", markersize=5, linestyle="none")
          axis.plot(1/R1_mesh_data[i], 1/yield_mesh_data_means[i], "-o", color=ammonium_colors[i])

     # Plot yield 1D scans over ammonium
     axis = figure.add_subplot(1, 2, 2)
     axis.text(-0.3, 1.05, "B", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Reciprocal added ammonium\nconcentration (1/mM)", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("Reciprocal growth yield (1/OD 600 nm)", fontsize=colimitation_plots.axis_label_size)
     axis.set_xlim([-0.2, 6.2])
     axis.set_xticks(numpy.linspace(0, 6, 7))
     axis.set_ylim([0, 80])

     norm = colors.SymLogNorm(linthresh=R1_mesh_data.T[:, 0][2], vmin=min(R1_mesh_data.T[:, 0]), vmax=max(R1_mesh_data.T[:, 0]))
     cmap = matplotlib.colormaps["cividis"]
     sm = pyplot.cm.ScalarMappable(cmap=cmap, norm=norm)
     sm.set_array([])
     glucose_colors = sm.to_rgba(R1_mesh_data.T[:, 0])
     colorbar = figure.colorbar(sm, ax=axis)
     colorbar.set_label("Added glucose\nconcentration (mM)", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=25)
     for i in range(len(yield_mesh_data_means.T)):
          #axis.errorbar(R2_mesh_data.T[i], yield_mesh_data_means.T[i], yerr=yield_mesh_data_sds.T[i], color=glucose_colors[i], marker="o", markersize=5, linestyle="none")
          axis.plot(1/R2_mesh_data.T[i], 1/yield_mesh_data_means.T[i], "-o", color=glucose_colors[i])

     figure.savefig("Figure_S22.pdf", bbox_inches="tight")

################################################################################


if __name__ == '__main__':
     main()
