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

     # Get parameters for a rate model fitted to the experimental data
     rate_model_to_plot = "pat_rmin"
     dataframe = pandas.read_csv("./Data/Rate_data_fits_rep123.csv", index_col=False)
     dataframe_filtered = dataframe[dataframe["Model name"] == rate_model_to_plot]
     row_index = dataframe_filtered.index[0]
     gmax, a1, a2, R1min, R2min = [float(x) for x in str(dataframe_filtered["Parameter values"][row_index]).split(";")]

     # Set up meshes for plots of fitted model
     R1_list_rate = numpy.linspace(1e-4, 0.1, 100)
     R2_list_rate = numpy.linspace(1e-4, 0.025, 100)
     R1_mesh_rate, R2_mesh_rate = numpy.meshgrid(R1_list_rate, R2_list_rate)
     L1_mesh_rate = numpy.zeros(R1_mesh_rate.shape)
     L2_mesh_rate = numpy.zeros(R1_mesh_rate.shape)
     Meff_mesh_rate = numpy.zeros(R1_mesh_rate.shape)

     # Iterate over points in R1-R2 space and calculate limitation coefficients and Meff
     for i in range(len(R2_list_rate)):
          for j in range(len(R1_list_rate)):
               L1_mesh_rate[i][j] = colimitation_models.CalcLimCoeff(0, [R1_list_rate[j], R2_list_rate[i]], [a1, a2], gmax, rate_model_to_plot.strip("_rmin"))
               L1_mesh_rate[i][j] = colimitation_models.CalcLimCoeff(1, [R1_list_rate[j], R2_list_rate[i]], [a1, a2], gmax, rate_model_to_plot.strip("_rmin"))
               Meff_mesh_rate[i][j] = colimitation_models.CalcMeff([R1_list_rate[j], R2_list_rate[i]], [a1, a2], gmax, rate_model_to_plot.strip("_rmin"))

################################################################################

     # Get parameters for a yield model fitted to the experimental data
     yield_model_to_plot = "pat"
     dataframe = pandas.read_csv("./Data/Yield_data_fits_rep123.csv", index_col=False)
     dataframe_filtered = dataframe[dataframe["Model name"] == yield_model_to_plot]
     row_index = dataframe_filtered.index[0]
     Nmax, s1, s2 = [float(x) for x in str(dataframe_filtered["Parameter values"][row_index]).split(";")]

     # Set up meshes for plots of fitted model
     R1_list_yield = numpy.linspace(1e-4, 10, 100)
     R2_list_yield = numpy.linspace(1e-4, 10, 100)
     R1_mesh_yield, R2_mesh_yield = numpy.meshgrid(R1_list_yield, R2_list_yield)
     L1_mesh_yield = numpy.zeros(R1_mesh_yield.shape)
     L2_mesh_yield = numpy.zeros(R1_mesh_yield.shape)
     Meff_mesh_yield = numpy.zeros(R1_mesh_yield.shape)

     # Iterate over points in R1-R2 space and calculate limitation coefficients and Meff
     for i in range(len(R2_list_yield)):
          for j in range(len(R1_list_yield)):
               L1_mesh_yield[i][j] = colimitation_models.CalcLimCoeff(0, [R1_list_yield[j], R2_list_yield[i]], [s1, s2], Nmax, yield_model_to_plot.strip("_rmin"))
               L1_mesh_yield[i][j] = colimitation_models.CalcLimCoeff(1, [R1_list_yield[j], R2_list_yield[i]], [s1, s2], Nmax, yield_model_to_plot.strip("_rmin"))
               Meff_mesh_yield[i][j] = colimitation_models.CalcMeff([R1_list_yield[j], R2_list_yield[i]], [s1, s2], Nmax, yield_model_to_plot.strip("_rmin"))

################################################################################

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Initialize figure
     figure = pyplot.figure(figsize=(4*2, 3*1))
     figure.subplots_adjust(wspace=0.6, hspace=0.4)

################################################################################

     # Number of rate-limiting resources
     axis = figure.add_subplot(1, 2, 1)
     axis.text(-0.33, 1.05, "A", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Glucose concentration (mM)", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("Ammonium concentration (mM)", fontsize=colimitation_plots.axis_label_size)

     contourfs = axis.contourf(R1_mesh_rate, R2_mesh_rate, Meff_mesh_rate, cmap="plasma", levels=numpy.linspace(1, 3, 11), zorder=0)
     contours = axis.contour(R1_mesh_rate, R2_mesh_rate, Meff_mesh_rate, levels=contourfs.levels, colors="black", linewidths=0.25, zorder=1)
     colorbar = figure.colorbar(contourfs, ax=axis, shrink=1.0)
     colorbar.set_label("Effective number of rate-\nlimiting factors $M^\mathrm{rate}_\mathrm{eff}$", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=25)

     axis.set_xlim([0, max(R1_list_rate)])
     axis.set_ylim([0, max(R2_list_rate)])
     axis.set_xticks(numpy.linspace(*axis.get_xlim(), 6))
     axis.set_yticks(numpy.linspace(*axis.get_ylim(), 6))

################################################################################

     # Number of yield-limiting resources
     axis = figure.add_subplot(1, 2, 2)
     axis.text(-0.25, 1.05, "B", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Glucose concentration (mM)", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("Ammonium concentration (mM)", fontsize=colimitation_plots.axis_label_size)

     contourfs = axis.contourf(R1_mesh_yield, R2_mesh_yield, Meff_mesh_yield, cmap="plasma", levels=numpy.linspace(1, 3, 11), zorder=0)
     contours = axis.contour(R1_mesh_yield, R2_mesh_yield, Meff_mesh_yield, levels=contourfs.levels, colors="black", linewidths=0.25, zorder=1)
     colorbar = figure.colorbar(contourfs, ax=axis, shrink=1.0)
     colorbar.set_label("Effective number of yield-\nlimiting factors $M^\mathrm{yield}_\mathrm{eff}$", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=25)

     axis.set_xlim([0, max(R1_list_yield)])
     axis.set_ylim([0, max(R2_list_yield)])
     axis.set_xticks(numpy.linspace(*axis.get_xlim(), 6))
     axis.set_yticks(numpy.linspace(*axis.get_ylim(), 6))

################################################################################

     figure.savefig("Figure_S31.pdf", bbox_inches="tight")

################################################################################


if __name__ == '__main__':
     main()
