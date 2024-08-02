import numpy
import matplotlib
from matplotlib import pyplot
import colimitation_plots
import colimitation_models


def main():

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Set stoichiometry
     Y1 = 1.5
     Y2 = 1

     # Initialize figure     
     figure = pyplot.figure(figsize=(4*3, 3*1))
     figure.subplots_adjust(wspace=0.3, hspace=0.2)

################################################################################

     # Plot schematic phase portrait for constant stoichiometry where growth 
     # stops at zero resources

     axis = figure.add_subplot(1, 3, 1)
     axis.text(-0.1, 1.05, "A", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Resource concentration $R_1$")
     axis.set_ylabel("Resource concentration $R_2$")
     axis.set_xlim([0, 1])
     axis.set_ylim([0, 1])
     pyplot.setp(axis.get_xticklabels(), visible=False)
     pyplot.setp(axis.get_yticklabels(), visible=False)
     axis.set_xticks([])
     axis.set_yticks([])

     R1_range = numpy.linspace(1, 0, 100)
     R2_range = (Y1/Y2)*R1_range
     for intercept in numpy.arange(1, -2, -0.2):
          line, = axis.plot(R1_range, R2_range + intercept, color="0.8")
          colimitation_plots.add_arrow_to_line2D(axis, line, arrow_locs=numpy.linspace(0, 1, 10), arrowstyle='-|>')
     line, = axis.plot(R1_range, R2_range, color="black")
     colimitation_plots.add_arrow_to_line2D(axis, line, arrow_locs=numpy.linspace(0, 1, 10), arrowstyle='-|>')
     axis.text(0.27, 0.8, r"\textbf{Trajectories that}" + "\n" + r"\textbf{exhaust resource 1}" + "\n" + r"\textbf{first} ($\Omega_1$)", transform=axis.transAxes, ha="center", va="center", fontsize=colimitation_plots.axis_label_size)
     axis.text(0.6, 0.25, r"\textbf{Trajectories that exhaust}" + "\n" + r"\textbf{resource 2 first} ($\Omega_2$)", transform=axis.transAxes, ha="center", va="center", fontsize=colimitation_plots.axis_label_size)

################################################################################

     # Plot schematic phase portrait for constant stoichiometry where growth 
     # stops at nonzero resources

     axis = figure.add_subplot(1, 3, 2)
     axis.text(-0.1, 1.05, "B", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Resource concentration $R_1$")
     axis.set_ylabel("Resource concentration $R_2$")
     axis.set_xlim([0, 1])
     axis.set_ylim([0, 1])
     pyplot.setp(axis.get_xticklabels(), visible=False)
     pyplot.setp(axis.get_yticklabels(), visible=False)
     axis.set_xticks([])
     axis.set_yticks([])

     # Plot ZNGI when growth stops at nonzero resource (e.g., when birth rate = death rate), 
     # assuming additive growth rate model
     gmax = 1
     d = 0.1
     R1_range_contour = numpy.linspace(1/(gmax/d - 1), 1, 100)
     R2_range_contour = [R1/((gmax/d - 1)*R1 - 1)  for R1 in R1_range_contour]
     axis.plot(R1_range_contour, R2_range_contour, "--", color="gray")
     axis.text(0.15, 0.9, "$g = 0$ at\n" + r"$R_i \neq 0$", transform=axis.transAxes, color="black", ha="left", va="top", fontsize=colimitation_plots.axis_label_size)

     # Plot trajectories ending on ZNGI
     R1initial_range = numpy.linspace(1/(gmax/d - 1), 2, 20)
     R2initial_range = 2 - R1initial_range
     K1, K2 = 1, 1
     for i in range(len(R1initial_range)):
          R1initial = R1initial_range[i]
          R2initial = R2initial_range[i]
          R1star, R2star = colimitation_models.CalcRStarAdditiveModel(gmax, K1, K2, d, R1initial, R2initial, Y1, Y2)
          line, = axis.plot(numpy.linspace(R1initial, R1star, 100), numpy.linspace(R2initial, R2star, 100), color="0.8")
          colimitation_plots.add_arrow_to_line2D(axis, line, arrow_locs=numpy.arange(0, 1, 0.2), arrowstyle='-|>')

################################################################################

     # Plot schematic phase portrait for variable stoichiometry where growth 
     # stops at zero resources

     axis = figure.add_subplot(1, 3, 3)
     axis.text(-0.1, 1.05, "C", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Resource concentration $R_1$")
     axis.set_ylabel("Resource concentration $R_2$")
     axis.set_xlim([0, 1])
     axis.set_ylim([0, 1])
     pyplot.setp(axis.get_xticklabels(), visible=False)
     pyplot.setp(axis.get_yticklabels(), visible=False)
     axis.set_xticks([])
     axis.set_yticks([])

     R1_range = numpy.linspace(1, 0, 100)
     R2_range = numpy.sqrt(R1_range)
     for intercept in numpy.arange(1, -2, -0.2):
          line, = axis.plot(R1_range, R2_range + intercept, color="0.8")
          colimitation_plots.add_arrow_to_line2D(axis, line, arrow_locs=numpy.linspace(0, 1, 10), arrowstyle='-|>')
     line, = axis.plot(R1_range, R2_range, color="black")
     colimitation_plots.add_arrow_to_line2D(axis, line, arrow_locs=numpy.linspace(0, 1, 10), arrowstyle='-|>')
     axis.text(0.27, 0.8, r"\textbf{Trajectories that}" + "\n" + r"\textbf{exhaust resource 1}" + "\n" + r"\textbf{first} ($\Omega_1$)", transform=axis.transAxes, ha="center", va="center", fontsize=colimitation_plots.axis_label_size)
     axis.text(0.6, 0.25, r"\textbf{Trajectories that exhaust}" + "\n" + r"\textbf{resource 2 first} ($\Omega_2$)", transform=axis.transAxes, ha="center", va="center", fontsize=colimitation_plots.axis_label_size)

################################################################################

     figure.savefig("Figure_S4.pdf", bbox_inches="tight")

################################################################################


if __name__ == '__main__':
     main()
