import numpy
import matplotlib
from matplotlib import pyplot
import colimitation_models
import consumer_resource_model
import colimitation_plots


def main():

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Set default growth rate model
     growth_rate_model = "pat"

     # Initialize figure     
     figure = pyplot.figure(figsize=(3.5*3, 3*2))
     figure.subplots_adjust(hspace=0.4, wspace=0.7)

################################################################################

     # Resource, trait, limitation, and colimitation under Blackman model
     Rs_blackman= [0, 0.99, 1.01, 5]
     zs_blackman = [0, 0.99, 1, 1]
     Ls_blackman = [1, 1, 0, 0]
     Meff_blackman = [1, 1, 1, 1]

     # Resource, trait, limitation, and colimitation under Monod model
     Rs_monod = numpy.linspace(0, 5, 50)
     zs_monod = Rs_monod/(Rs_monod + 1)
     Ls_monod = 1/(Rs_monod + 1)
     Meff_monod = 1/numpy.maximum(Ls_monod, 1 - Ls_monod)

     # Plot trait dependence on resource concentration
     axis = pyplot.subplot2grid((14, 3), (0, 0), rowspan=2)
     axis.set_ylabel("Growth\ntrait\n(rate or\nyield)", fontsize=colimitation_plots.axis_label_size)
     axis.set_xlim([0, max(Rs_blackman)])
     axis.set_ylim([-0.1, 1.1])
     axis.set_xticks([])
     pyplot.setp(axis.get_xticklabels(), visible=False)
     axis.set_yticks([])
     pyplot.setp(axis.get_yticklabels(), visible=False)
     axis.plot(Rs_blackman, zs_blackman, "-", color="tab:cyan")
     axis.plot(Rs_monod, zs_monod, "--", color="tab:cyan")

     # Plot limitation coefficient dependence on resource concentration
     axis = pyplot.subplot2grid((14, 3), (2, 0), rowspan=2)
     axis.set_ylabel("Limitation\ncoefficient\n$L$", fontsize=colimitation_plots.axis_label_size)
     axis.set_xlim([0, max(Rs_blackman)])
     axis.set_ylim([-0.1, 1.1])
     axis.set_xticks([])
     pyplot.setp(axis.get_xticklabels(), visible=False)
     axis.tick_params(axis="y", labelsize=0.8*colimitation_plots.tick_label_size)
     axis.plot(Rs_blackman, Ls_blackman, "-", color="tab:red")
     axis.plot(Rs_monod, Ls_monod, "--", color="tab:red")

     # Plot colimitation dependence on resource concentration
     axis = pyplot.subplot2grid((14, 3), (4, 0), rowspan=2)
     axis.set_xlabel("Resource concentration $R$", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("No.\nlimiting\nresources\n$M_\mathrm{eff}$", fontsize=colimitation_plots.axis_label_size)
     axis.set_xlim([0, max(Rs_blackman)])
     axis.set_ylim([0.9, 2.1])
     axis.set_xticks([])
     pyplot.setp(axis.get_xticklabels(), visible=False)
     axis.tick_params(axis="y", labelsize=0.8*colimitation_plots.tick_label_size)
     axis.plot(Rs_blackman, Meff_blackman, "-", color="black")
     axis.plot(Rs_monod, Meff_monod, "--", color="black")

################################################################################

     # Initialize axis for limitation state phase diagram
     axis = pyplot.subplot2grid((14, 3), (0, 1), rowspan=6)
     axis.text(-2, 1.08, "A", transform=axis.transAxes, fontsize=1.3*colimitation_plots.panel_label_size)
     axis.text(-0.25, 1.08, "B", transform=axis.transAxes, fontsize=1.3*colimitation_plots.panel_label_size)
     axis.text(1.4, 1.08, "C", transform=axis.transAxes, fontsize=1.3*colimitation_plots.panel_label_size)
     axis.set_xlabel("Resource concentration $R_1$", fontsize=1.2*colimitation_plots.axis_label_size)
     axis.set_ylabel("Resource concentration $R_2$", fontsize=1.2*colimitation_plots.axis_label_size)
     axis.set_xlim([0, 3])
     axis.set_ylim([0, 3])
     axis.set_xticks([])
     axis.set_yticks([])
     pyplot.setp(axis.get_xticklabels(), visible=False)
     pyplot.setp(axis.get_yticklabels(), visible=False)
     axis.set_aspect("equal")

     # Region for R1 single limitation
     cmap = matplotlib.colormaps["Blues"]
     axis.fill_betweenx([0, 1, 3], [0, 1, 1], x2=0, color=cmap(0.5), alpha=0.5)
     axis.text(0.45, 2, r"\textbf{Single}" + "\n" + r"\textbf{lim}" + "\n" + "$R_1$", fontsize=1.3*colimitation_plots.axis_label_size, ha="center", va="center")

     # Region for R2 single limitation
     cmap = matplotlib.colormaps["Oranges"]
     axis.fill_between([0, 1, 3], [0, 1, 1], y2=0, color=cmap(0.5), alpha=0.5)
     axis.text(2, 0.45, r"\textbf{Single lim} $R_2$", fontsize=1.3*colimitation_plots.axis_label_size, ha="center", va="center")

     # Region for R3 single limitation
     cmap = matplotlib.colormaps["Greens"]
     axis.fill_between([1, 3], [1, 1], y2=3, color=cmap(0.5), alpha=0.5)
     axis.text(2.25, 2.5, r"\textbf{Single lim}" + "\n" + r"\textbf{implicit}" + "\n" + r"\textbf{factors}", fontsize=1.3*colimitation_plots.axis_label_size, ha="center", va="center")

     # Region for R1,R3 colimitation
     axis.fill_betweenx([1, 3], [0.9, 0.9], x2=[1.1, 1.1], color="yellow")
     axis.text(1, 2, r"\textbf{Colim} $R_1$,\textbf{imp}", rotation=90, fontsize=0.9*colimitation_plots.axis_label_size, ha="center", va="center")

     # Region for R2,R3 colimitation
     axis.fill_between([1, 3], [0.9, 0.9], y2=[1.1, 1.1], color="yellow")
     axis.text(2, 1, r"\textbf{Colim} $R_2$,\textbf{imp}", fontsize=0.9*colimitation_plots.axis_label_size, ha="center", va="center")

     # Region for R1,R2 colimitation
     axis.fill_between([0, 1.1], [0 - 0.2/numpy.sqrt(2), 1.1 - 0.2/numpy.sqrt(2)], y2=[0 + 0.2/numpy.sqrt(2), 1.1 + 0.2/numpy.sqrt(2)], color="yellow")
     axis.text(0.5, 0.5, r"\textbf{Colim} $R_1$,$R_2$", rotation=45, fontsize=0.9*colimitation_plots.axis_label_size, ha="center", va="center")

     # Region for R1,R2,R3 colimitation
     axis.annotate(r"\textbf{Colim}" + "\n" + r"$R_1$,$R_2$,\textbf{imp}", xytext=(1.75, 1.5), xy=(1, 1), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.tick_label_size, color="black"), color="black", ha="center", va="center", size=0.9*colimitation_plots.axis_label_size)

################################################################################

     # Calculate over batch growth cycle: biomass and resource concentrations, limitation coefficients, and colimitation

     # Initialize consumer-resource model
     gmax = 1
     gmaxes = numpy.array([gmax])
     Ks = numpy.array([[1, 2]])
     Ys = numpy.array([1e2, 2e2])
     production_consumption = numpy.array([-1/Ys])
     R0s = numpy.array([1, 1])
     N0s = numpy.array([1])
     myModel = consumer_resource_model.ConsumerResourceModel(gmaxes, Ks, production_consumption, [growth_rate_model])

     # Set saturation time
     tsat = 40

     # Time points during growth cycle at which to calculate concentrations
     time_points_growth = numpy.linspace(0, tsat, 200)

     # Calculate biomass and resource concentrations over one batch growth cycle
     Ns, Rs = myModel.CalcBatchDynamics(N0s, R0s, time_points_growth)

     # Calculate instantaneous limitation coefficient trajectories (shape is resources, time points)
     L1s = [colimitation_models.CalcLimCoeff(0, Rs.T[t], gmax/Ks[0], gmax, growth_rate_model) for t in range(len(time_points_growth))]
     L2s = [colimitation_models.CalcLimCoeff(1, Rs.T[t], gmax/Ks[0], gmax, growth_rate_model) for t in range(len(time_points_growth))]

     # Calculate instantaneous Meff over time points
     Meffs = numpy.array([colimitation_models.CalcMeff(Rs.T[t], gmax/Ks[0], gmax, growth_rate_model) for t in range(len(time_points_growth))])

################################################################################

     # Plot biomass concentration over batch cycle on left axis
     axis = pyplot.subplot2grid((14, 3), (0, 2), rowspan=3)
     axis.set_xticklabels([])
     axis.set_ylabel("Biomass\nconcentration $N$", fontsize=colimitation_plots.axis_label_size)
     axis.set_yscale("log")
     axis.set_xlim([min(time_points_growth), max(time_points_growth)])
     axis.set_ylim([10**(-0.1), 10**(2.1)])
     axis.tick_params(axis="y", labelsize=0.8*colimitation_plots.tick_label_size, pad=0) 
     axis.plot(time_points_growth, Ns[0], "-", color="black", linewidth=2)
     axis.text(9, 7, "$N$", fontsize=colimitation_plots.axis_label_size, color="black")
     axis.set_title(r"\textbf{Batch growth}", fontsize=1.3*colimitation_plots.axis_label_size)

     # Plot resource concentrations over batch cycle on right axis
     axis_twinx = axis.twinx()
     axis_twinx.set_ylabel("Resource\nconcentration $R_i$", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=20)
     axis_twinx.set_yscale("log")
     axis_twinx.set_ylim([10**(-4.2), 10**(0.2)])
     axis_twinx.set_yticks(numpy.logspace(-4, 0, 5))
     axis_twinx.tick_params(axis="y", labelsize=0.8*colimitation_plots.tick_label_size, pad=0) 
     axis_twinx.plot(time_points_growth, Rs[0], "-", color="blue", linewidth=2)
     axis_twinx.plot(time_points_growth, Rs[1], "-", color="orange", linewidth=2)
     axis_twinx.text(23, 1e-3, "$R_1$", fontsize=colimitation_plots.axis_label_size, color="blue")
     axis_twinx.text(32, 1e-1, "$R_2$", fontsize=colimitation_plots.axis_label_size, color="orange")

     # Plot number of limiting factors over batch cycle on left axis
     axis = pyplot.subplot2grid((14, 3), (3, 2), rowspan=3)
     axis.set_xlabel("Time", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("No. rate-limiting\nfactors $M^\mathrm{rate}_\mathrm{eff}$", fontsize=colimitation_plots.axis_label_size)
     axis.set_xlim([min(time_points_growth), max(time_points_growth)])
     axis.set_ylim([0.9, 3.1])
     axis.set_yticks(numpy.linspace(1, 3, 5))
     axis.tick_params(axis="y", labelsize=0.8*colimitation_plots.tick_label_size, pad=0) 
     axis.plot(time_points_growth, Meffs, "-", color="black", linewidth=2)
     axis.text(10, 2.55, "$M^\mathrm{rate}_\mathrm{eff}$", fontsize=colimitation_plots.axis_label_size, color="black")

     # Plot limitation coefficients for all resources over batch cycle on right axis
     axis_twinx = axis.twinx()
     axis_twinx.set_ylabel("Rate limitation\ncoefficient $L^\mathrm{rate}_i$", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=22)
     axis_twinx.set_ylim([-0.05, 1.05])
     axis_twinx.set_yticks(numpy.linspace(0, 1, 6))
     axis_twinx.tick_params(axis="y", labelsize=0.8*colimitation_plots.tick_label_size, pad=0) 
     axis_twinx.plot(time_points_growth, L1s, "-", color="blue", linewidth=2)
     axis_twinx.plot(time_points_growth, L2s, "-", color="orange", linewidth=2)
     axis_twinx.text(30, 0.83, "$L^\mathrm{rate}_1$", fontsize=colimitation_plots.axis_label_size, color="blue")
     axis_twinx.text(3, 0.6, "$L^\mathrm{rate}_2$", fontsize=colimitation_plots.axis_label_size, color="orange")

################################################################################

     # Initialize axis for resource deplation trajectory over batch cycle
     axis = pyplot.subplot2grid((14, 3), (8, 1), rowspan=6)
     axis.text(-2, 1.08, "D", transform=axis.transAxes, fontsize=1.3*colimitation_plots.panel_label_size)
     axis.text(-0.25, 1.08, "E", transform=axis.transAxes, fontsize=1.3*colimitation_plots.panel_label_size)
     axis.text(1.4, 1.08, "F", transform=axis.transAxes, fontsize=1.3*colimitation_plots.panel_label_size)
     axis.set_xlabel("Resource concentration $R_1$")
     axis.set_ylabel("Resource concentration $R_2$")
     axis.set_xlim([0, 1.5])
     axis.set_ylim([0, 1.5])
     axis.set_xticks(numpy.linspace(0, 1.5, 4))
     axis.set_yticks(numpy.linspace(0, 1.5, 4))
     axis.set_aspect("equal")
     axis.set_title(r"\textbf{Batch growth}", fontsize=1.3*colimitation_plots.axis_label_size)

     # Draw resource depletion trajectory
     line, = axis.plot(Rs[0], Rs[1], "-", linewidth=3, color="tab:purple", zorder=0)
     colimitation_plots.add_arrow_to_line2D(axis, line, arrow_locs=[0.3, 0.55, 0.8], arrowsize=1.5)

     # Mark initial and final points
     axis.scatter(R0s[0], R0s[1], s=40, c="black", zorder=1)
     axis.annotate("Initial concentrations\n$R_1(0)$, $R_2(0)$", xytext=(0.95, 1.25), xy=(R0s[0], R0s[1]), arrowprops=dict(arrowstyle="-|>", color="black"), ha="center", fontsize=colimitation_plots.axis_label_size)
     axis.scatter(Rs[0][-1], Rs[1][-1], s=40, marker="s", c="black", zorder=1)
     axis.annotate("Final concentrations\n($R_1$ is depleted)", xytext=(0.7, 0.1), xy=(Rs[0][-1], Rs[1][-1]), arrowprops=dict(arrowstyle="-|>", color="black"), ha="center", fontsize=colimitation_plots.axis_label_size)

     # Annotate slope
     axis.text(0.45, 0.65, "Slope =\n$R_2$:$R_1$ stoichiometry", rotation=numpy.arctan(0.5)*180/numpy.pi, ha="center", fontsize=colimitation_plots.axis_label_size)

################################################################################

     # Calculate for chemostat: biomass and resource concentrations, limitation coefficients, and colimitation

     # Initialize consumer-resource model for chemostat: same as batch above but
     # with constant resource supply concentrations and dilution rate
     Rsources = numpy.array([1, 1])
     dilution_rate = 0.2

     # Initial biomass and resource concentrations
     R0s = Rsources
     N0s = numpy.array([1])

     # Set saturation time
     tsat = 200

     # Time points at which to calculate concentrations
     time_points_growth = numpy.linspace(0, tsat, 200)

     # Calculate biomass and resource concentrations for chemostat
     Ns, Rs = myModel.CalcChemostatDynamics(N0s, R0s, time_points_growth, Rsources, dilution_rate)

     # Calculate instantaneous limitation coefficient trajectories (shape is resources, time points)
     L1s = [colimitation_models.CalcLimCoeff(0, Rs.T[t], gmax/Ks[0], gmax, growth_rate_model) for t in range(len(time_points_growth))]
     L2s = [colimitation_models.CalcLimCoeff(1, Rs.T[t], gmax/Ks[0], gmax, growth_rate_model) for t in range(len(time_points_growth))]

     # Calculate instantaneous Meff over time points
     Meffs = numpy.array([colimitation_models.CalcMeff(Rs.T[t], gmax/Ks[0], gmax, growth_rate_model) for t in range(len(time_points_growth))])

################################################################################

     # Plot biomass concentration for chemostat on left axis
     axis = pyplot.subplot2grid((14, 3), (8, 0), rowspan=3)
     axis.set_xticklabels([])
     axis.set_ylabel("Biomass\nconcentration $N$", fontsize=colimitation_plots.axis_label_size)
     axis.set_yscale("log")
     axis.set_xlim([min(time_points_growth), max(time_points_growth)])
     axis.set_ylim([10**(-0.1), 10**(2.1)])
     axis.tick_params(axis="y", labelsize=0.8*colimitation_plots.tick_label_size, pad=0) 
     axis.plot(time_points_growth, Ns[0], "-", color="black", linewidth=2)
     axis.text(20, 3, "$N$", fontsize=colimitation_plots.axis_label_size, color="black")
     axis.set_title(r"\textbf{Chemostat growth}", fontsize=1.3*colimitation_plots.axis_label_size)

     # Plot resource concentrations for chemostat on right axis
     axis_twinx = axis.twinx()
     axis_twinx.set_ylabel("Resource\nconcentration $R_i$", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=20)
     axis_twinx.set_ylim([0.37, 1.03])
     #axis_twinx.set_yticks([0.6, 0.7, 0.8, 0.9, 1.0])
     axis_twinx.tick_params(axis="y", labelsize=0.8*colimitation_plots.tick_label_size, pad=0) 
     axis_twinx.plot(time_points_growth, Rs[0], "-", color="blue", linewidth=2)
     axis_twinx.plot(time_points_growth, Rs[1], "-", color="orange", linewidth=2)
     axis_twinx.text(125, 0.45, "$R_1$", fontsize=colimitation_plots.axis_label_size, color="blue")
     axis_twinx.text(175, 0.75, "$R_2$", fontsize=colimitation_plots.axis_label_size, color="orange")

     # Plot number of limiting factors in chemostat on left axis
     axis = pyplot.subplot2grid((14, 3), (11, 0), rowspan=3)
     axis.set_xlabel("Time", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("No. rate-limiting\nfactors $M^\mathrm{rate}_\mathrm{eff}$", fontsize=colimitation_plots.axis_label_size)
     axis.set_xlim([min(time_points_growth), max(time_points_growth)])
     axis.set_ylim([0.9, 3.1])
     axis.set_yticks(numpy.linspace(1, 3, 5))
     axis.tick_params(axis="y", labelsize=0.8*colimitation_plots.tick_label_size, pad=0) 
     axis.plot(time_points_growth, Meffs, "-", color="black", linewidth=2)
     axis.text(100, 2.4, "$M^\mathrm{rate}_\mathrm{eff}$", fontsize=colimitation_plots.axis_label_size, color="black")

     # Plot limitation coefficients for all resources in chemostat on right axis
     axis_twinx = axis.twinx()
     axis_twinx.set_ylabel("Rate limitation\ncoefficient $L^\mathrm{rate}_i$", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=25)
     axis_twinx.set_ylim([-0.05, 1.05])
     axis_twinx.set_yticks(numpy.linspace(0, 1, 6))
     axis_twinx.tick_params(axis="y", labelsize=0.8*colimitation_plots.tick_label_size, pad=0) 
     axis_twinx.plot(time_points_growth, L1s, "-", color="blue", linewidth=2)
     axis_twinx.plot(time_points_growth, L2s, "-", color="orange", linewidth=2)
     axis_twinx.text(100, 0.18, "$L^\mathrm{rate}_1$", fontsize=colimitation_plots.axis_label_size, color="blue")
     axis_twinx.text(10, 0.6, "$L^\mathrm{rate}_2$", fontsize=colimitation_plots.axis_label_size, color="orange")

################################################################################

     # Initialize axis for resource deplation trajectory over batch cycle
     axis = pyplot.subplot2grid((14, 3), (8, 2), rowspan=6)
     axis.set_xlabel("Resource concentration $R_1$")
     axis.set_ylabel("Resource concentration $R_2$")
     axis.set_xlim([0, 1.5])
     axis.set_ylim([0, 1.5])
     axis.set_xticks(numpy.linspace(0, 1.5, 4))
     axis.set_yticks(numpy.linspace(0, 1.5, 4))
     axis.set_aspect("equal")
     axis.set_title(r"\textbf{Chemostat growth}", fontsize=1.3*colimitation_plots.axis_label_size)

     # Draw resource depletion trajectory
     line, = axis.plot(Rs[0], Rs[1], "-", linewidth=3, color="tab:purple", zorder=0)
     colimitation_plots.add_arrow_to_line2D(axis, line, arrow_locs=[1/3 + 0.1, 2/3 + 0.1], arrowsize=1.5)

     # Mark initial and final points
     axis.scatter(Rsources[0], Rsources[1], s=40, c="black", zorder=1)
     axis.annotate("Source concentrations\n$R_1^\mathrm{source}$, $R_2^\mathrm{source}$", xytext=(0.95, 1.25), xy=(Rsources[0], Rsources[1]), arrowprops=dict(arrowstyle="-|>", color="black"), ha="center", fontsize=colimitation_plots.axis_label_size)
     axis.scatter(Rs[0][-1], Rs[1][-1], s=40, marker="s", c="black", zorder=1)
     axis.annotate("Zero net-growth isocline\n(ZNGI): $g(R_1,R_2) = d$", xytext=(0.7, 0.1), xy=(Rs[0][-1], Rs[1][-1]), arrowprops=dict(arrowstyle="-|>", color="black"), ha="center", fontsize=colimitation_plots.axis_label_size)

     # Calculate and draw ZNGI
     R1_list_model = numpy.linspace(0, 1.5, 100)
     R2_list_model = numpy.linspace(0, 1.5, 100)
     R1_mesh_model, R2_mesh_model = numpy.meshgrid(R1_list_model, R2_list_model)
     rate_mesh_model = colimitation_models.CalcTraitForFit((R1_mesh_model, R2_mesh_model), [gmax, gmax/Ks[0][0], gmax/Ks[0][1]], growth_rate_model)
     axis.contour(R1_mesh_model, R2_mesh_model, rate_mesh_model, levels=[dilution_rate, 1], linestyles=["--"], colors="gray", zorder=-1)

################################################################################

     figure.savefig("Figure_1.pdf", bbox_inches="tight")

################################################################################


if __name__ == '__main__':
     main()
