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

     # List of models for which to plot fits
     models_to_plot = ["liebig_blackman", "liebig_monod", "pat"]

################################################################################

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Initialize figure
     figure = pyplot.figure(figsize=(4*len(models_to_plot), 3*2))
     figure.subplots_adjust(wspace=0.6, hspace=0.5)

################################################################################

     for m in range(len(models_to_plot)):

          # Get parameters for a model fitted to the experimental data
          model_to_plot = models_to_plot[m]
          dataframe = pandas.read_csv("./Data/Yield_data_fits_rep123.csv", index_col=False)
          dataframe_filtered = dataframe[dataframe["Model name"] == model_to_plot]
          row_index = dataframe_filtered.index[0]
          params = [float(x) for x in str(dataframe_filtered["Parameter values"][row_index]).split(";")]
     
          # Plot yield 1D scans over glucose
          axis = figure.add_subplot(2, len(models_to_plot), m + 1)
          axis.text(-0.3, 1.05, chr(65 + m), transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
          axis.set_xlabel("Added glucose concentration (mM)", fontsize=colimitation_plots.axis_label_size)
          axis.set_ylabel("Growth yield (OD 600 nm)", fontsize=colimitation_plots.axis_label_size)
          axis.set_ylim([-0.01, 0.6])
          axis.set_title(r"\textbf{" + colimitation_models.all_2D_trait_models_formatted[colimitation_models.all_2D_trait_models.index(model_to_plot)] + "}", fontsize=1.5*colimitation_plots.axis_label_size, y = 1.2)
     
          norm = colors.SymLogNorm(linthresh=R2_mesh_data[:, 0][1], vmin=min(R2_mesh_data[:, 0]), vmax=max(R2_mesh_data[:, 0]))
          cmap = matplotlib.colormaps["viridis"]
          sm = pyplot.cm.ScalarMappable(cmap=cmap, norm=norm)
          sm.set_array([])
          ammonium_colors = sm.to_rgba(R2_mesh_data[:, 0])
          colorbar = figure.colorbar(sm, ax=axis)
          colorbar.set_label("Added ammonium\nconcentration (mM)", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=25)
          for i in range(len(R1_mesh_data)):
               for r in range(num_replicates):
                    xys = numpy.array([[R1_mesh_data[i][j], yield_mesh_data[i, j, r]] for j in range(len(R1_mesh_data[i])) if not isinstance(yield_mesh_data[i, j, r], str)])
                    if len(xys) != 0:
                         axis.scatter(xys.T[0], xys.T[1], color=ammonium_colors[i], marker="o", s=5)
               R1s_model_fit = numpy.concatenate((numpy.linspace(R1_mesh_data[i][0], R1_mesh_data[i][1], 20), numpy.logspace(numpy.log10(R1_mesh_data[i][1]), numpy.log10(max(R1_mesh_data[i])), 50)))
               yields_model_fit = colimitation_models.CalcTraitForFit((R1s_model_fit, R2_mesh_data[i][0]), params, model_to_plot)
               axis.plot(R1s_model_fit, yields_model_fit, "-", color=ammonium_colors[i])
          axis.set_xscale("symlog", linthresh=R1_mesh_data[0][2])

          # Plot yield 1D scans over ammonium
          axis = figure.add_subplot(2, len(models_to_plot), m + 1 + len(models_to_plot))
          axis.text(-0.3, 1.05, chr(65 + m + len(models_to_plot)), transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
          axis.set_xlabel("Added ammonium concentration (mM)", fontsize=colimitation_plots.axis_label_size)
          axis.set_ylabel("Growth yield (OD 600 nm)", fontsize=colimitation_plots.axis_label_size)
          axis.set_ylim([-0.01, 0.6])
     
          norm = colors.SymLogNorm(linthresh=R1_mesh_data.T[:, 0][2], vmin=min(R1_mesh_data.T[:, 0]), vmax=max(R1_mesh_data.T[:, 0]))
          cmap = matplotlib.colormaps["cividis"]
          sm = pyplot.cm.ScalarMappable(cmap=cmap, norm=norm)
          sm.set_array([])
          glucose_colors = sm.to_rgba(R1_mesh_data.T[:, 0])
          colorbar = figure.colorbar(sm, ax=axis)
          colorbar.set_label("Added glucose\nconcentration (mM)", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=25)
          for i in range(len(R2_mesh_data.T)):
               for r in range(num_replicates):
                    xys = numpy.array([[R2_mesh_data[j][i], yield_mesh_data[j, i, r]] for j in range(len(R2_mesh_data[:, i])) if not isinstance(yield_mesh_data[j, i, r], str)])
                    if len(xys) != 0:
                         axis.scatter(xys.T[0], xys.T[1], color=glucose_colors[i], marker="o", s=5)
               R2s_model_fit = numpy.concatenate((numpy.linspace(R2_mesh_data.T[i][0], R2_mesh_data.T[i][1], 20), numpy.logspace(numpy.log10(R2_mesh_data.T[i][1]), numpy.log10(max(R2_mesh_data.T[i])), 50)))
               yields_model_fit = colimitation_models.CalcTraitForFit((R1_mesh_data.T[i][0], R2s_model_fit), params, model_to_plot)
               axis.plot(R2s_model_fit, yields_model_fit, "-", color=glucose_colors[i])
          axis.set_xscale("symlog", linthresh=R2_mesh_data.T[0][1])

################################################################################

     figure.savefig("Figure_S22.pdf", bbox_inches="tight")

################################################################################


if __name__ == '__main__':
     main()
