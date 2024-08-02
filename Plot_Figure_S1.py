import numpy
from matplotlib import pyplot
import matplotlib
import colimitation_models
import colimitation_plots


def main():

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Initialize figure
     figure = pyplot.figure(figsize=(4*4, 2.5))
     figure.subplots_adjust(hspace=0.1, wspace=0.2)

     # Set up parameters for scan
     ss = [1, 1]
     gmax = 1
     R1_range = numpy.linspace(1e-3, 10, 100)
     R2_range = numpy.linspace(1e-3, 10, 100)
     R1_mesh, R2_mesh = numpy.meshgrid(R1_range, R2_range)

     # Set model name
     model = "pat"

     # Initialize growth mesh
     g_mesh = numpy.zeros((len(R2_range), len(R1_range)))

     # Iterate over points in R1-R2 space
     for i in range(len(R2_range)):
          for j in range(len(R1_range)):

               # Set R1, R2 point
               R1, R2 = R1_range[j], R2_range[i]
               Rs = [R1, R2]

               # Calculate growth 
               g_mesh[i][j] = colimitation_models.CalcTraitForFit(Rs, [gmax, *ss], model)

################################################################################

     # 2D scan with supplementations on different backgrounds

     # Set up axis
     axis = pyplot.subplot2grid((2, 8), (0, 0), colspan=2, rowspan=2)
     axis.text(-0.1, 1.01, "A", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.text(1.3, 1.01, "B", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlim([0, max(R1_range)])
     axis.set_ylim([0, max(R2_range)])
     pyplot.setp(axis.get_xticklabels(), visible=False)
     pyplot.setp(axis.get_yticklabels(), visible=False)
     axis.set_xticks([])
     axis.set_yticks([])
     axis.set_xlabel("Resource 1 concentration", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("Resource 2 concentration", fontsize=colimitation_plots.axis_label_size)

     # Plot contours
     contourfs = axis.contourf(R1_mesh, R2_mesh, g_mesh/(1 - g_mesh), cmap="Oranges", alpha=0.75) #, levels=growth_levels)
     contours = axis.contour(R1_mesh, R2_mesh, g_mesh/(1 - g_mesh), levels=contourfs.levels, colors="black", linewidths=0.25)
     colorbar = figure.colorbar(contourfs, ax=axis, shrink=1.0)
     colorbar.set_label("Growth trait (rate or yield)", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=15)
     colorbar.set_ticks([])
     axis.set_aspect("equal")

     supplementation_scenarios = [(8, 0.2, 2, 2), (3, 0.5, 2, 2), (6.5, 6.5, 2, 2), (0.2, 0.2, 2, 2)]
     scenario_labels = [ r"\textbf{I: Single}" + "\n" + r"\textbf{lim}", 
                         r"\textbf{II: Serial}" + "\n" + r"\textbf{lim}", 
                         r"\textbf{III: Add}" + "\n" + r"\textbf{colim}", 
                         r"\textbf{IV: Super}" + "\n" + r"\textbf{colim}"]
     scenario_titles = [r"\textbf{I}", r"\textbf{II}", r"\textbf{III}", r"\textbf{IV}"]

     for s in range(len(supplementation_scenarios)):

          R1, R2, dR1, dR2 = supplementation_scenarios[s]
          g0 = colimitation_models.CalcTraitForFit([R1, R2], [gmax, *ss], model)
          g1 = colimitation_models.CalcTraitForFit([R1 + dR1, R2], [gmax, *ss], model)
          g2 = colimitation_models.CalcTraitForFit([R1, R2 + dR2], [gmax, *ss], model)
          g12 = colimitation_models.CalcTraitForFit([R1 + dR1, R2 + dR2], [gmax, *ss], model)
          transform = lambda x: (x - g0)/(g12 - g0) + 0.1  

          # Plot supplementation on main axis
          axis.scatter(R1, R2, color="black", zorder=10)
          axis.annotate("", xytext=(R1, R2), xy=(R1 + dR1, R2), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.axis_label_size, color="tab:blue", linewidth=2))
          axis.annotate("", xytext=(R1, R2), xy=(R1, R2 + dR2), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.axis_label_size, color="tab:red", linewidth=2))
          axis.annotate("", xytext=(R1, R2), xy=(R1 + dR1, R2 + dR2), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.axis_label_size, color="tab:purple", linewidth=2))
          axis.text(R1 + 0.5*dR1, R2 + dR2, scenario_titles[s], fontsize=colimitation_plots.axis_label_size, ha="center", va="bottom")

          # Set up sub axis to plot bars
          axis_sub = pyplot.subplot2grid((2, 8), (int(s/2), 2 + s%2), colspan=1, rowspan=1)
          axis_sub.set_ylim([0, 1.2])
          axis_sub.set_xlim([-0.5, 3.5])
          axis_sub.set_xticks([0, 1, 2, 3])
          if s == 0 or s == 1:
               axis_sub.set_xticklabels([])
          else:
               axis_sub.set_xticklabels(["No\nsupp", "+R1", "+R2", "+R1\n+R2"])
          axis_sub.set_yticks([])
          if s == 0 or s == 2:
               axis_sub.set_ylabel("Growth trait", fontsize=colimitation_plots.axis_label_size)
          axis_sub.text(0.02, 0.95, scenario_labels[s], transform=axis_sub.transAxes, fontsize=colimitation_plots.axis_label_size, ha="left", va="top")
          axis_sub.bar([0, 1, 2, 3], transform(numpy.array([g0, g1, g2, g12])), color=["black", "tab:blue", "tab:red", "tab:purple"])

################################################################################

     # 2D scan with different supplementations on the same background: 
     # additive colimitation can look like single limitation

     # Set up axis
     axis = pyplot.subplot2grid((2, 8), (0, 4), colspan=2, rowspan=2)
     axis.text(-0.1, 1.01, "C", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlim([0, max(R1_range)])
     axis.set_ylim([0, max(R2_range)])
     pyplot.setp(axis.get_xticklabels(), visible=False)
     pyplot.setp(axis.get_yticklabels(), visible=False)
     axis.set_xticks([])
     axis.set_yticks([])
     axis.set_xlabel("Resource 1 concentration", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("Resource 2 concentration", fontsize=colimitation_plots.axis_label_size)

     # Plot contours
     contourfs = axis.contourf(R1_mesh, R2_mesh, g_mesh/(1 - g_mesh), cmap="Oranges", alpha=0.75) #, levels=growth_levels)
     contours = axis.contour(R1_mesh, R2_mesh, g_mesh/(1 - g_mesh), levels=contourfs.levels, colors="black", linewidths=0.25)
     colorbar = figure.colorbar(contourfs, ax=axis, shrink=1.0)
     colorbar.set_label("Growth trait (rate or yield)", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=15)
     colorbar.set_ticks([])
     axis.set_aspect("equal")

     # Get additive scenario but modify dR1, dR2
     R1, R2, dR1, dR2 = supplementation_scenarios[2]
     dR1, dR2 = 0.01, 2.5
     dR_offset = 1
     g0 = colimitation_models.CalcTraitForFit([R1, R2], [gmax, *ss], model)
     g1 = colimitation_models.CalcTraitForFit([R1 + dR1, R2], [gmax, *ss], model)
     g2 = colimitation_models.CalcTraitForFit([R1, R2 + dR2], [gmax, *ss], model)
     g12 = colimitation_models.CalcTraitForFit([R1 + dR1, R2 + dR2], [gmax, *ss], model)
     transform = lambda x: (x - g0)/(g12 - g0) + 0.1  

     # Plot supplementation on main axis
     axis.scatter(R1, R2, color="black", zorder=10)
     axis.annotate("", xytext=(R1, R2), xy=(R1 + dR1 + dR_offset, R2), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.axis_label_size, color="tab:blue", linewidth=2))
     axis.annotate("", xytext=(R1, R2), xy=(R1, R2 + dR2 + dR_offset), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.axis_label_size, color="tab:red", linewidth=2))
     axis.annotate("", xytext=(R1, R2), xy=(R1 + dR1 + dR_offset, R2 + dR2 + dR_offset), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.axis_label_size, color="tab:purple", linewidth=2))
     axis.text(R1 - 0.2, R2 + 0.5*(dR2 + dR_offset), r"\textbf{V}", fontsize=colimitation_plots.axis_label_size, ha="right", va="center")

     # Plot inset
     axis_inset = axis.inset_axes([0.1, 0.18, 0.7, 0.4])
     axis_inset.set_ylim([0, 1.2])
     axis_inset.set_xlim([-0.5, 3.5])
     axis_inset.set_xticks([0, 1, 2, 3])
     axis_inset.set_xticklabels([r"\textbf{No}" + "\n" + r"\textbf{supp}", r"\textbf{+R1}", r"\textbf{+R2}", r"\textbf{+R1}" + "\n" + r"\textbf{+R2}"], fontsize=colimitation_plots.legend_label_size)
     axis_inset.set_yticks([])
     axis_inset.set_ylabel(r"\textbf{Growth rate}", fontsize=colimitation_plots.legend_label_size)
     axis_inset.text(0.02, 0.95, r"\textbf{V: Single}" + "\n" + r"\textbf{lim}", transform=axis_inset.transAxes, fontsize=colimitation_plots.axis_label_size, ha="left", va="top")
     axis_inset.bar([0, 1, 2, 3], transform(numpy.array([g0, g1, g2, g12])), color=["black", "tab:blue", "tab:red", "tab:purple"])

################################################################################

     # 2D scan with different supplementations on the same background: 
     # Super additive colimitation can look like additive colimitation

     # Set up axis
     axis = pyplot.subplot2grid((2, 8), (0, 6), colspan=2, rowspan=2)
     axis.text(-0.1, 1.01, "D", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlim([0, max(R1_range)])
     axis.set_ylim([0, max(R2_range)])
     pyplot.setp(axis.get_xticklabels(), visible=False)
     pyplot.setp(axis.get_yticklabels(), visible=False)
     axis.set_xticks([])
     axis.set_yticks([])
     axis.set_xlabel("Resource 1 concentration", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("Resource 2 concentration", fontsize=colimitation_plots.axis_label_size)

     # Plot contours
     contourfs = axis.contourf(R1_mesh, R2_mesh, g_mesh/(1 - g_mesh), cmap="Oranges", alpha=0.75) #, levels=growth_levels)
     contours = axis.contour(R1_mesh, R2_mesh, g_mesh/(1 - g_mesh), levels=contourfs.levels, colors="black", linewidths=0.25)
     colorbar = figure.colorbar(contourfs, ax=axis, shrink=1.0)
     colorbar.set_label("Growth rate (rate or yield)", fontsize=colimitation_plots.axis_label_size, rotation=270, labelpad=15)
     colorbar.set_ticks([])
     axis.set_aspect("equal")

     # Get additive scenario but modify dR1, dR2
     R1, R2, dR1, dR2 = supplementation_scenarios[3]
     dR1, dR2 = 0.01, 0.01
     dR_offset = 1
     g0 = colimitation_models.CalcTraitForFit([R1, R2], [gmax, *ss], model)
     g1 = colimitation_models.CalcTraitForFit([R1 + dR1, R2], [gmax, *ss], model)
     g2 = colimitation_models.CalcTraitForFit([R1, R2 + dR2], [gmax, *ss], model)
     g12 = colimitation_models.CalcTraitForFit([R1 + dR1, R2 + dR2], [gmax, *ss], model)
     transform = lambda x: (x - g0)/(g12 - g0) + 0.1  

     # Plot supplementation on main axis
     axis.scatter(R1, R2, color="black", zorder=10)
     axis.annotate("", xytext=(R1, R2), xy=(R1 + dR1 + dR_offset, R2), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.axis_label_size, color="tab:blue", linewidth=2))
     axis.annotate("", xytext=(R1, R2), xy=(R1, R2 + dR2 + dR_offset), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.axis_label_size, color="tab:red", linewidth=2))
     axis.annotate("", xytext=(R1, R2), xy=(R1 + dR1 + dR_offset, R2 + dR2 + dR_offset), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.axis_label_size, color="tab:purple", linewidth=2))
     axis.text(R1 + 0.5*(dR1 + dR_offset), R2 + dR2 + dR_offset, r"\textbf{VI}", fontsize=colimitation_plots.axis_label_size, ha="center", va="bottom")

     # Plot inset
     axis_inset = axis.inset_axes([0.25, 0.55, 0.7, 0.4])
     axis_inset.set_ylim([0, 1.2])
     axis_inset.set_xlim([-0.5, 3.5])
     axis_inset.set_xticks([0, 1, 2, 3])
     axis_inset.set_xticklabels([r"\textbf{No}" + "\n" + r"\textbf{supp}", r"\textbf{+R1}", r"\textbf{+R2}", r"\textbf{+R1}" + "\n" + r"\textbf{+R2}"], fontsize=colimitation_plots.legend_label_size)
     axis_inset.set_yticks([])
     axis_inset.set_ylabel(r"\textbf{Growth trait}", fontsize=colimitation_plots.legend_label_size)
     axis_inset.text(0.02, 0.95, r"\textbf{VI: Add}" + "\n" + r"\textbf{colim}", transform=axis_inset.transAxes, fontsize=colimitation_plots.axis_label_size, ha="left", va="top")
     axis_inset.bar([0, 1, 2, 3], transform(numpy.array([g0, g1, g2, g12])), color=["black", "tab:blue", "tab:red", "tab:purple"])

################################################################################

     figure.savefig("Figure_S1.pdf", bbox_inches="tight")

################################################################################


if __name__ == '__main__':
     main()
