import ternary
import matplotlib
from matplotlib import pyplot
import numpy
import colimitation_plots


def main():

     # Meff function
     Meff = lambda x: 1/max(x[0], x[1], x[2])
     
     # Resolution of grid for heatmap 
     scale = 20

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)
     
     # Initialize figure and axis
     figure = pyplot.figure(figsize=(5.5, 4))
     axis = figure.add_subplot(111)
     axis.axis("off")
     figure, ternary_axis = ternary.figure(ax=axis, scale=scale)
     
     # Plot heatmap
     ternary_axis.heatmapf(Meff, boundary=True, cmap=matplotlib.colormaps["plasma"], cbarlabel="Effective number of limiting factors $M_\mathrm{eff}$", vmin=1.0, vmax=3.0)    

     # Draw boundaries
     ternary_axis.boundary(linewidth=2.0, axes_colors={"b": "black", "l": "black", "r": "black"})
     
     # Write axis labels
     ternary_axis.left_axis_label("Limitation coefficient $L_1$", offset=0.16, color="black")
     ternary_axis.right_axis_label("Limitation coefficient $L_2$", offset=0.16, color="black")
     ternary_axis.bottom_axis_label("Limitation coefficient $L_3$", offset=0.06, color="black")
     
     # Set and format axes ticks
     ternary_axis.ticks(ticks=[i/10 for i in range(10 + 1)], axis='rlb', linewidth=1, clockwise=True, axes_colors={"b": "black", "l": "black", "r": "black"}, offset=0.02, tick_formats="%0.1f")
     
     # Save figure
     ternary_axis.clear_matplotlib_ticks()
     ternary_axis._redraw_labels()
     figure.savefig("Figure_S2.pdf", bbox_inches="tight")

################################################################################


if __name__ == '__main__':
     main()