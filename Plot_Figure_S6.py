import numpy
import matplotlib
from matplotlib import pyplot
import colimitation_plots
import colimitation_models


def main():

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Initialize figure     
     figure = pyplot.figure(figsize=(4*2, 3*3))
     figure.subplots_adjust(wspace=0.3, hspace=0.5)

################################################################################

     # Plot model with resource consumption for maintenance 

     # Define parameters for this model
     gmax = 2
     s1main = 0.5
     s2main = 10
     s1growth = 2
     s2growth = 0.1
     a1, a2 = 2, 2
     R1max = 1

     # Define stoichiometry functions for this model
     s1_function = lambda R1, R2, gmax, a1, a2, s1main, s2main, s1growth, s2growth: (1/s1growth + 1/s1main/colimitation_models.CalcTraitForFit([R1, R2], [gmax, a1, a2], "additive"))**(-1)
     s2_function = lambda R1, R2, gmax, a1, a2, s1main, s2main, s1growth, s2growth: (1/s2growth + 1/s2main/colimitation_models.CalcTraitForFit([R1, R2], [gmax, a1, a2], "additive"))**(-1)
     dR2dR1 = lambda R1, R2, gmax, a1, a2, s1main, s2main, s1growth, s2growth: s1_function(R1, R2, gmax, a1, a2, s1main, s2main, s1growth, s2growth)/s2_function(R1, R2, gmax, a1, a2, s1main, s2main, s1growth, s2growth)

     # Calculate phase line separating Omega_1 and Omega_2 (regions where 
     # resource 1 or resource 2 is first exhausted
     phase_line = colimitation_models.CalcResourceDepletionPhaseLine(dR2dR1, [gmax, a1, a2, s1main, s2main, s1growth, s2growth], R1max)

     # Initialize meshes for resource concentrations, stoichiometry, and total yield
     R1_range = numpy.linspace(1e-3, 1, 100)
     R2_range = numpy.linspace(1e-3, 1, 100)
     R1_mesh, R2_mesh = numpy.meshgrid(R1_range, R2_range)
     stoich1_mesh = numpy.zeros((len(R2_range), len(R1_range)))
     stoich2_mesh = numpy.zeros((len(R2_range), len(R1_range)))
     N_mesh = numpy.zeros((len(R2_range), len(R1_range)))

     # Iterate over points in R1-R2 space
     for i in range(len(R2_range)):
          for j in range(len(R1_range)):

               # Set R1, R2 point
               R1, R2 = R1_range[j], R2_range[i]

               # Calculate stoichiometry
               stoich1_mesh[i][j] = -1/s1_function(R1, R2, gmax, a1, a2, s1main, s2main, s1growth, s2growth)
               stoich2_mesh[i][j] = -1/s2_function(R1, R2, gmax, a1, a2, s1main, s2main, s1growth, s2growth)

               # Calculate total biomass yield
               N_mesh[i][j] = colimitation_models.CalcYieldVariableStoichiometry(R1, R2, s1_function, s2_function, phase_line, (gmax, a1, a2, s1main, s2main, s1growth, s2growth))

     # Plot phase portrait
     axis = figure.add_subplot(3, 2, 1)
     axis.text(-0.3, 1.05, "A", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Resource concentration $R_1$")
     axis.set_ylabel("Resource concentration $R_2$")
     axis.set_xlim([0, 1])
     axis.set_ylim([0, 1])
     axis.set_aspect("equal")
     axis.text(-1, 0.5, r"\textbf{Batch dynamics}" + "\n" + r"\textbf{with maintenance}" + "\n" + r"\textbf{resource}" + "\n" + r"\textbf{consumption}", fontsize=1.5*colimitation_plots.axis_label_size, transform=axis.transAxes, ha="center", va="center")

     # Plot stream lines for resource depletion
     axis.streamplot(R1_mesh, R2_mesh, stoich1_mesh, stoich2_mesh, color="0.7", linewidth=0.5, arrowsize=0.5, density=1)
     axis.text(0.35, 0.75, r"\textbf{Trajectories that}" + "\n" + r"\textbf{exhaust resource 1}" + "\n" + r"\textbf{first} ($\Omega_1$)", transform=axis.transAxes, ha="center", va="center", fontsize=colimitation_plots.legend_label_size)
     axis.text(0.78, 0.2, r"\textbf{Trajectories}" + "\n" + r"\textbf{that exhaust}" + "\n" + r"\textbf{resource 2}" + "\n" + r"\textbf{first} ($\Omega_2$)", transform=axis.transAxes, ha="center", va="center", fontsize=colimitation_plots.legend_label_size)

     # Plot phase line
     line, = axis.plot(R1_range[::-1], phase_line(R1_range[::-1]), "-", color="black")
     colimitation_plots.add_arrow_to_line2D(axis, line, arrow_locs=numpy.arange(0, 1, 0.2), arrowstyle='-|>')

     # Plot total biomass yield
     axis = figure.add_subplot(3, 2, 2)
     axis.text(-0.3, 1.05, "B", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Resource concentration $R_1$")
     axis.set_ylabel("Resource concentration $R_2$")
     axis.set_xlim([0, 1])
     axis.set_ylim([0, 1])
     axis.set_aspect("equal")
     contourfs = axis.contourf(R1_mesh, R2_mesh, N_mesh, cmap="viridis") #, levels=numpy.linspace(0, 1, 11))
     colorbar = figure.colorbar(contourfs, ax=axis, shrink=1.0)
     colorbar.set_label("Growth yield $N$", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=15)
     line, = axis.plot(R1_range[::-1], phase_line(R1_range[::-1]), "-", color="white")
     colimitation_plots.add_arrow_to_line2D(axis, line, arrow_locs=numpy.arange(0, 1, 0.2), arrowstyle='-|>')
     axis.text(0.35, 0.75, r"\textbf{Trajectories that}" + "\n" + r"\textbf{exhaust resource 1}" + "\n" + r"\textbf{first} ($\Omega_1$)", transform=axis.transAxes, ha="center", va="center", fontsize=colimitation_plots.legend_label_size, color="white")
     axis.text(0.78, 0.2, r"\textbf{Trajectories}" + "\n" + r"\textbf{that exhaust}" + "\n" + r"\textbf{resource 2}" + "\n" + r"\textbf{first} ($\Omega_2$)", transform=axis.transAxes, ha="center", va="center", fontsize=colimitation_plots.legend_label_size, color="white")

################################################################################

     # Plot model with resource consumption that varies with proteome

     # Define parameters for this model
     gmax = 1
     theta_min = 0
     theta_max = 1
     s1A = 1
     s1B = 1
     s2A = 10
     s2B = 0.1
     a1, a2 = 1, 1
     R1max = 50

     # Define stoichiometry functions for this model
     theta = lambda R1, R2, gmax, a1, a2, theta_min, theta_max: colimitation_models.CalcTraitForFit([R1, R2], [gmax, a1, a2], "additive")*(theta_max - theta_min)/gmax + theta_min
     s1_function = lambda R1, R2, gmax, K1, K2, theta_min, theta_max, s1A, s1B, s2A, s2B: (theta(R1, R2, gmax, a1, a2, theta_min, theta_max)/s1A + (1 - theta(R1, R2, gmax, a1, a2, theta_min, theta_max))/s1B)**(-1)
     s2_function = lambda R1, R2, gmax, K1, K2, theta_min, theta_max, s1A, s1B, s2A, s2B: (theta(R1, R2, gmax, a1, a2, theta_min, theta_max)/s2A + (1 - theta(R1, R2, gmax, a1, a2, theta_min, theta_max))/s2B)**(-1)
     dR2dR1 = lambda R1, R2, gmax, a1, a2, theta_min, theta_max, s1A, s1B, s2A, s2B: s1_function(R1, R2, gmax, a1, a2, theta_min, theta_max, s1A, s1B, s2A, s2B)/s2_function(R1, R2, gmax, a1, a2, theta_min, theta_max, s2A, s2B, s2A, s2B)

     # Calculate phase line separating Omega_1 and Omega_2 (regions where 
     # resource 1 or resource 2 is first exhausted
     phase_line = colimitation_models.CalcResourceDepletionPhaseLine(dR2dR1, [gmax, a1, a2, theta_min, theta_max, s1A, s1B, s2A, s2B], R1max)

     # Initialize meshes for resource concentrations, stoichiometry, and total yield
     R1_range = numpy.linspace(1e-3, 50, 100)
     R2_range = numpy.linspace(1e-3, 50, 100)
     R1_mesh, R2_mesh = numpy.meshgrid(R1_range, R2_range)
     stoich1_mesh = numpy.zeros((len(R2_range), len(R1_range)))
     stoich2_mesh = numpy.zeros((len(R2_range), len(R1_range)))
     N_mesh = numpy.zeros((len(R2_range), len(R1_range)))

     # Iterate over points in R1-R2 space
     for i in range(len(R2_range)):
          for j in range(len(R1_range)):

               # Set R1, R2 point
               R1, R2 = R1_range[j], R2_range[i]

               # Calculate stoichiometry
               stoich1_mesh[i][j] = -1/s1_function(R1, R2, gmax, a1, a2, theta_min, theta_max, s1A, s1B, s2A, s2B)
               stoich2_mesh[i][j] = -1/s2_function(R1, R2, gmax, a1, a2, theta_min, theta_max, s1A, s1B, s2A, s2B)

               # Calculate total biomass yield
               N_mesh[i][j] = colimitation_models.CalcYieldVariableStoichiometry(R1, R2, s1_function, s2_function, phase_line, (gmax, a1, a2, theta_min, theta_max, s1A, s1B, s2A, s2B))

     # Plot phase portrait
     axis = figure.add_subplot(3, 2, 3)
     axis.text(-0.3, 1.05, "C", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Resource concentration $R_1$")
     axis.set_ylabel("Resource concentration $R_2$")
     axis.set_xlim([0, 50])
     axis.set_ylim([0, 50])
     axis.set_aspect("equal")
     axis.text(-1, 0.5, r"\textbf{Batch dynamics}" + "\n" + r"\textbf{with dynamic}" + "\n" + r"\textbf{proteome}" + "\n" + r"\textbf{allocation}", fontsize=1.5*colimitation_plots.axis_label_size, transform=axis.transAxes, ha="center", va="center")
     axis.set_xticks(numpy.linspace(0, 50, 6))
     axis.set_yticks(numpy.linspace(0, 50, 6))

     # Plot stream lines for resource depletion
     axis.streamplot(R1_mesh, R2_mesh, stoich1_mesh, stoich2_mesh, color="0.7", linewidth=0.5, arrowsize=0.5, density=1)
     axis.text(0.22, 0.99, r"\textbf{Trajectories}" + "\n" + r"\textbf{that exhaust}" + "\n" + r"\textbf{resource 1}" + "\n" + r"\textbf{first} ($\Omega_1$)", transform=axis.transAxes, ha="center", va="top", fontsize=colimitation_plots.legend_label_size)
     axis.text(0.6, 0.4, r"\textbf{Trajectories}" + "\n" + r"\textbf{that exhaust}" + "\n" + r"\textbf{resource 2}" + "\n" + r"\textbf{first} ($\Omega_2$)", transform=axis.transAxes, ha="center", va="center", fontsize=colimitation_plots.legend_label_size)

     # Plot phase line
     line, = axis.plot(R1_range[::-1], phase_line(R1_range[::-1]), "-", color="black")
     colimitation_plots.add_arrow_to_line2D(axis, line, arrow_locs=numpy.arange(0, 1, 0.2), arrowstyle='-|>')

     # Plot total biomass yield
     axis = figure.add_subplot(3, 2, 4)
     axis.text(-0.3, 1.05, "D", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Resource concentration $R_1$")
     axis.set_ylabel("Resource concentration $R_2$")
     axis.set_xlim([0, 50])
     axis.set_ylim([0, 50])
     axis.set_aspect("equal")
     contourfs = axis.contourf(R1_mesh, R2_mesh, N_mesh, cmap="viridis") #, levels=numpy.linspace(0, 1, 11))
     colorbar = figure.colorbar(contourfs, ax=axis, shrink=1.0)
     colorbar.set_label("Growth yield $N$", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=15)
     line, = axis.plot(R1_range[::-1], phase_line(R1_range[::-1]), "-", color="white")
     colimitation_plots.add_arrow_to_line2D(axis, line, arrow_locs=numpy.arange(0, 1, 0.2), arrowstyle='-|>')
     axis.text(0.22, 0.99, r"\textbf{Trajectories}" + "\n" + r"\textbf{that exhaust}" + "\n" + r"\textbf{resource 1}" + "\n" + r"\textbf{first} ($\Omega_1$)", transform=axis.transAxes, ha="center", va="top", fontsize=colimitation_plots.legend_label_size, color="white")
     axis.text(0.6, 0.4, r"\textbf{Trajectories}" + "\n" + r"\textbf{that exhaust}" + "\n" + r"\textbf{resource 2}" + "\n" + r"\textbf{first} ($\Omega_2$)", transform=axis.transAxes, ha="center", va="center", fontsize=colimitation_plots.legend_label_size, color="white")
     axis.set_xticks(numpy.linspace(0, 50, 6))
     axis.set_yticks(numpy.linspace(0, 50, 6))

################################################################################

     # Plot model where growth stops at nonzero concentrations

     # Define parameters for this model
     gmax = 1
     a1, a2 = 1, 1
     d = 0.1 
     s1, s2 = 1.5, 1

     # Initialize meshes for resource concentrations, stoichiometry, and total yield
     R1_range = numpy.linspace(1e-3, 1, 100)
     R2_range = numpy.linspace(1e-3, 1, 100)
     R1_mesh, R2_mesh = numpy.meshgrid(R1_range, R2_range)
     stoich1_mesh = numpy.zeros((len(R2_range), len(R1_range)))
     stoich2_mesh = numpy.zeros((len(R2_range), len(R1_range)))
     N_mesh = numpy.zeros((len(R2_range), len(R1_range)))

     # Iterate over points in R1-R2 space
     for i in range(len(R2_range)):
          for j in range(len(R1_range)):

               # Set R1, R2 point
               R1, R2 = R1_range[j], R2_range[i]

               # Calculate stoichiometry
               stoich1_mesh[i][j] = -1/s1
               stoich2_mesh[i][j] = -1/s2              

               # Calculate total biomass yield
               N_mesh[i][j] = colimitation_models.CalcYieldChemostatAdditiveModel(gmax, a1, a2, d, R1, R2, s1, s2)

     # Plot phase portrait
     axis = figure.add_subplot(3, 2, 5)
     axis.text(-0.3, 1.05, "E", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Resource concentration $R_1$")
     axis.set_ylabel("Resource concentration $R_2$")
     axis.set_xlim([0, 1])
     axis.set_ylim([0, 1])
     axis.set_aspect("equal")
     axis.text(-1, 0.5, r"\textbf{Chemostat dynamics}" + "\n" + r"\textbf{(growth stops}" + "\n" + r"\textbf{at nonzero}" + "\n" + r"\textbf{resource}" + "\n" + r"\textbf{concentration)}", fontsize=1.5*colimitation_plots.axis_label_size, transform=axis.transAxes, ha="center", va="center")

     # Plot stream lines for resource depletion
     streams = axis.streamplot(R1_mesh, R2_mesh, stoich1_mesh, stoich2_mesh, color="0.7", linewidth=0.5, arrowsize=0.5, density=1)

     # Plot ZNGI (mask stream lines for growth rates below the lower limit)
     R1min_ZNGI = 1/(gmax/d - 1)
     R1_range_ZNGI= numpy.linspace(R1min_ZNGI, 1, 100)
     R2_range_ZNGI = [R1/(R1/R1min_ZNGI - 1) for R1 in R1_range_ZNGI]
     axis.fill_between([0, R1min_ZNGI], [0, 0], y2=1, color="white", zorder=2) 
     axis.fill_between(R1_range_ZNGI, R2_range_ZNGI, y2=0, color="white", zorder=2)
     axis.plot(R1_range_ZNGI, R2_range_ZNGI, "--", color="gray")
     axis.text(0.7, 0.1, "$g = 0$", transform=axis.transAxes, color="black", ha="left", va="top", fontsize=colimitation_plots.axis_label_size)

     # Plot total biomass yield
     axis = figure.add_subplot(3, 2, 6)
     axis.text(-0.3, 1.05, "F", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Resource concentration $R_1$")
     axis.set_ylabel("Resource concentration $R_2$")
     axis.set_xlim([0, 1])
     axis.set_ylim([0, 1])
     axis.set_aspect("equal")
     contourfs = axis.contourf(R1_mesh, R2_mesh, N_mesh, cmap="viridis", levels=numpy.linspace(0, 1, 11))
     colorbar = figure.colorbar(contourfs, ax=axis, shrink=1.0)
     colorbar.set_label("Growth yield $N$", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=15)
     axis.plot(R1_range_ZNGI, R2_range_ZNGI, "--", color="gray")

################################################################################

     figure.savefig("Figure_S6.pdf", bbox_inches="tight")

################################################################################


if __name__ == '__main__':
     main()
