import numpy
from matplotlib import pyplot
import matplotlib
import colimitation_models
import colimitation_plots


def main():

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Initialize fixed parameter values
     a1, a2 = 1, 1
     gmax = 1

     # Initialize meshes of R1 and R2
     R1_range_large = numpy.linspace(1e-3, 10, 100)
     R2_range_large = numpy.linspace(1e-3, 10, 100)
     R1_mesh_large, R2_mesh_large = numpy.meshgrid(R1_range_large, R2_range_large)
     R1_range_small = numpy.linspace(1e-3, 5, 100)
     R2_range_small = numpy.linspace(1e-3, 5, 100)
     R1_mesh_small, R2_mesh_small = numpy.meshgrid(R1_range_small, R2_range_small)

     # Growth rate models to test
     growth_rate_models = ["saito_sub", "mean_monod_sub", "chem_depend"]
     model_labels = [r"\textbf{Saito et al.}" + "\n" + r"\textbf{substitutable model}", r"\textbf{Mean Monod}" + "\n" + r"\textbf{substitutable model}", r"\textbf{Chemically-dependent}" + "\n" + r"\textbf{model}"]

     # Initialize contour levels for growth rate
     growth_rate_levels = [numpy.linspace(0, 20, 11)/(numpy.linspace(0, 20, 11) + 1)] + 2*[numpy.linspace(0, 10, 11)/(numpy.linspace(0, 10, 11) + 1)]

     # Initial contour levels for limitation coefficients
     limitation_levels = [([0.04, 0.06, 0.1, 1], [0.04, 0.06, 0.1, 1]), 
                        ([0.08, 0.12, 0.18, 1], [0.08, 0.12, 0.18, 1]),
                        ([0.18, 0.25, 0.4, 1], [0.08, 0.2, 1])]
     cmap = matplotlib.colormaps["Blues"]
     limitation_color_levels_1 = [cmap(c) for c in numpy.linspace(0.4, 1, len(limitation_levels[0][0]))]
     cmap = matplotlib.colormaps["Oranges"]
     limitation_color_levels_2 = [cmap(c) for c in numpy.linspace(0.4, 1, len(limitation_levels[0][0]))]

     # Initialize figure
     figure = pyplot.figure(figsize=(4*len(growth_rate_models), 3*3))
     figure.subplots_adjust(wspace=0.5, hspace=0.5)

     # Iterate over growth rate models
     for m in range(len(growth_rate_models)):
          model = growth_rate_models[m]

          # Mesh of growth rate
          g_mesh = numpy.zeros((len(R2_range_large), len(R1_range_large)))

          # Meshes of limitation coefficients
          L1_mesh = numpy.zeros((len(R2_range_large), len(R1_range_large)))
          L2_mesh = numpy.zeros((len(R2_range_large), len(R1_range_large)))

          # Mesh of Meff
          Meff_mesh = numpy.zeros((len(R2_range_small), len(R1_range_small)))

          # Iterate over points in R1-R2 space
          for i in range(len(R2_range_large)):
               for j in range(len(R1_range_large)):

                    # Set R1, R2 point
                    R1, R2 = R1_range_large[j], R2_range_large[i]

                    # Calculate growth rate
                    g_mesh[i][j] = colimitation_models.CalcTraitForFit([R1, R2], [gmax, a1, a2], model)

                    # Calculate limitation coefficients
                    L1_mesh[i][j] = colimitation_models.CalcLimCoeff(0, [R1, R2], [a1, a2], gmax, model)
                    L2_mesh[i][j] = colimitation_models.CalcLimCoeff(1, [R1, R2], [a1, a2], gmax, model)

          # Iterate over points in R1-R2 space
          for i in range(len(R2_range_small)):
               for j in range(len(R1_range_small)):

                    # Set R1, R2 point
                    R1, R2 = R1_range_small[j], R2_range_small[i]

                    # Calculate Meff
                    Meff_mesh[i][j] = colimitation_models.CalcMeff([R1, R2], [R1, R2], gmax, model)

          # Initial growth rate landscape
          axis = figure.add_subplot(3, len(growth_rate_models), m + 1)
          contourfs = axis.contourf(R1_mesh_large, R2_mesh_large, g_mesh, cmap="viridis", levels=growth_rate_levels[m])
          contours = axis.contour(R1_mesh_large, R2_mesh_large, g_mesh, levels=contourfs.levels, colors="black", linewidths=0.25)
          colorbar = figure.colorbar(contourfs, ax=axis, shrink=1.0)
          colorbar.set_label("Growth rate $g(R_1, R_2)$", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=15)
          axis.set_aspect("equal")
          axis.set_xlabel("Resource concentration $R_1$", fontsize=colimitation_plots.axis_label_size)
          axis.set_ylabel("Resource concentration $R_2$", fontsize=colimitation_plots.axis_label_size)
          axis.text(0.5, 1.4, model_labels[m], fontsize=1.5*colimitation_plots.axis_label_size, fontweight="bold", ha="center", va="top", transform=axis.transAxes)
          axis.set_xlim([min(R1_range_large), max(R1_range_large)])
          axis.set_ylim([min(R2_range_large), max(R2_range_large)])
          axis.set_xticks(numpy.linspace(0, max(R1_range_large), 6))
          axis.set_yticks(numpy.linspace(0, max(R1_range_large), 6))

          # Contours of limitation coefficients for each resource, with invisible colorbar to match alignment of other panels
          axis = figure.add_subplot(3, len(growth_rate_models), m + 1 + len(growth_rate_models))
          contours = axis.contourf(R1_mesh_large, R2_mesh_large, L1_mesh, levels=limitation_levels[m][0], alpha=0.5, colors=limitation_color_levels_1) #cmap="Blues") #colors="blue", 
          contours = axis.contourf(R1_mesh_large, R2_mesh_large, L2_mesh, levels=limitation_levels[m][1], alpha=0.5, colors=limitation_color_levels_2) #cmap="Oranges") #colors="orange", 
          contours = axis.contour(R1_mesh_large, R2_mesh_large, g_mesh, colors="white", alpha=0)
          colorbar = figure.colorbar(contours, ax=axis, shrink=1.0)
          colorbar.set_ticks([])
          colorbar.outline.set_visible(False)
          axis.set_aspect("equal")
          axis.set_xlabel("Resource concentration $R_1$", fontsize=colimitation_plots.axis_label_size)
          axis.set_ylabel("Resource concentration $R_2$", fontsize=colimitation_plots.axis_label_size)
          axis.set_xlim([min(R1_range_large), max(R1_range_large)])
          axis.set_ylim([min(R2_range_large), max(R2_range_large)])
          axis.set_xticks(numpy.linspace(0, max(R1_range_large), 6))
          axis.set_yticks(numpy.linspace(0, max(R1_range_large), 6))
          if m == 0:
               axis.text(4, 0.6, "$L^\mathrm{rate}_1 > " + str(limitation_levels[m][0][2]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center")
               axis.text(7.5, 2, "$L^\mathrm{rate}_1 > " + str(limitation_levels[m][0][1]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center")
               axis.text(7.5, 4.5, "$L^\mathrm{rate}_1 > " + str(limitation_levels[m][0][0]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center")
               axis.text(0.6, 4, "$L^\mathrm{rate}_2 > " + str(limitation_levels[m][1][2]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center", rotation=90)
               axis.text(2, 7.5, "$L^\mathrm{rate}_2 > " + str(limitation_levels[m][1][1]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center", rotation=90)
               axis.text(4.5, 7.5, "$L^\mathrm{rate}_2 > " + str(limitation_levels[m][1][0]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center", rotation=90)
          elif m == 1:
               axis.text(4, 0.7, "$L^\mathrm{rate}_2 > " + str(limitation_levels[m][1][2]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center")
               axis.text(7.5, 2, "$L^\mathrm{rate}_2 > " + str(limitation_levels[m][1][1]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center")
               axis.text(7.5, 4.2, "$L^\mathrm{rate}_2 > " + str(limitation_levels[m][1][0]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center")
               axis.text(0.7, 4, "$L^\mathrm{rate}_1 > " + str(limitation_levels[m][0][2]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center", rotation=90)
               axis.text(2, 7.5, "$L^\mathrm{rate}_1 > " + str(limitation_levels[m][0][1]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center", rotation=90)
               axis.text(4.2, 7.5, "$L^\mathrm{rate}_1 > " + str(limitation_levels[m][0][0]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center", rotation=90)               
          else:
               axis.text(2, 0.5, "$L^\mathrm{rate}_2 > " + str(limitation_levels[m][1][1]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center")
               axis.text(4, 1.7, "$L^\mathrm{rate}_2 > " + str(limitation_levels[m][1][0]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center")
               axis.text(0.7, 7.5, "$L^\mathrm{rate}_1 > " + str(limitation_levels[m][0][2]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center", rotation=90)
               axis.text(2.4, 7.5, "$L^\mathrm{rate}_1 > " + str(limitation_levels[m][0][1]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center", rotation=90)
               axis.text(4.2, 7.5, "$L^\mathrm{rate}_1 > " + str(limitation_levels[m][0][0]) + "$", fontsize=colimitation_plots.axis_label_size, ha="center", va="center", rotation=90)               

          # Landscape of effective number of resources
          axis = figure.add_subplot(3, len(growth_rate_models), m + 1 + 2*len(growth_rate_models))
          contourfs = axis.contourf(R1_mesh_small, R2_mesh_small, Meff_mesh, cmap="plasma", levels=numpy.linspace(1, 3, 11))
          contours = axis.contour(R1_mesh_small, R2_mesh_small, Meff_mesh, levels=contourfs.levels, colors="black", linewidths=0.25)
          colorbar = figure.colorbar(contourfs, ax=axis, shrink=1.0)
          colorbar.set_label("Effective number of rate-\nlimiting factors $M^\mathrm{rate}_\mathrm{eff}$", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=25)
          axis.set_aspect("equal")
          axis.set_xlabel("Resource concentration $R_1$", fontsize=colimitation_plots.axis_label_size)
          axis.set_ylabel("Resource concentration $R_2$", fontsize=colimitation_plots.axis_label_size)
          axis.set_xlim([min(R1_range_small), max(R1_range_small)])
          axis.set_ylim([min(R2_range_small), max(R2_range_small)])
          axis.set_xticks(numpy.linspace(0, max(R1_range_small), 6))
          axis.set_yticks(numpy.linspace(0, max(R1_range_small), 6))

     figure.savefig("Figure_S4.pdf", bbox_inches="tight")

#################################################################################


if __name__ == '__main__':
     main()
