"""
This module contains parameters and a few functions to standardize plotting 
across figures.
"""


import numpy
import matplotlib.lines as mlines
import matplotlib.patches as mpatches


# Set font sizes
panel_label_size = 20
axis_label_size = 10
tick_label_size = 10
legend_label_size = 8
legend_label_size_small = 6


def SetrcParams(rcParams):
     """
     rcParams: Current state of rcParams from pyplot

     return: New state of rcParams with font settings
     """

     # Set fonts
     rcParams["font.family"] = "serif"
     rcParams["text.usetex"] = True
     rcParams["font.size"] = tick_label_size

     return rcParams


def CalcPointsLineSegment(xs, ys, num_segments):
     """
     xs: Pair of x coordinates
     ys: Pair of y coordinates
     num_segments: Number of segments to equally divide the whole line into

     return: List of (num_segments - 1) equally-spaced (x, y) points
     """

     # Unpack coordinates
     x1, x2 = xs
     y1, y2 = ys

     # Calculate slope
     slope = (y2 - y1)/(x2 - x1)

     # Calculate total distance along line segment
     Dtotal = numpy.sqrt((x2 - x1)**2 + (y2 - y1)**2)

     # Calculate equally-spaced points along line segment
     xpoints = [-Dtotal*(n/num_segments)/numpy.sqrt(1 + slope**2) + x1 for n in range(1, num_segments)]
     ypoints = [slope*(x - x1) + y1 for x in xpoints]

     # Combine x and y points into array of pairs
     xypoints = numpy.array([xpoints, ypoints])
     xypoints = xypoints.T

     return xypoints


def add_arrow_to_line2D(axes, line, arrow_locs=[0.2, 0.4, 0.6, 0.8], arrowstyle='-|>', arrowsize=1, transform=None):
     """
     Copied verbatim from https://stackoverflow.com/questions/26911898/matplotlib-curve-with-arrow-ticks
     
     Add arrows to a matplotlib.lines.Line2D at selected locations.
     
     Parameters:
     -----------
     axes: 
     line: Line2D object as returned by plot command
     arrow_locs: list of locations where to insert arrows, % of total length
     arrowstyle: style of the arrow
     arrowsize: size of the arrow
     transform: a matplotlib transform instance, default to data coordinates
     
     Returns:
     --------
     arrows: list of arrows
     """
     if not isinstance(line, mlines.Line2D):
          raise ValueError("expected a matplotlib.lines.Line2D object")
     x, y = line.get_xdata(), line.get_ydata()
     
     arrow_kw = {
        "arrowstyle": arrowstyle,
        "mutation_scale": 10 * arrowsize,
     }
     
     color = line.get_color()
     use_multicolor_lines = isinstance(color, numpy.ndarray)
     if use_multicolor_lines:
          raise NotImplementedError("multicolor lines not supported")
     else:
          arrow_kw['color'] = color
     
     linewidth = line.get_linewidth()
     if isinstance(linewidth, numpy.ndarray):
          raise NotImplementedError("multiwidth lines not supported")
     else:
          arrow_kw['linewidth'] = linewidth
     
     if transform is None:
          transform = axes.transData
     
     arrows = []
     for loc in arrow_locs:
          s = numpy.cumsum(numpy.sqrt(numpy.diff(x) ** 2 + numpy.diff(y) ** 2))
          n = numpy.searchsorted(s, s[-1] * loc)
          arrow_tail = (x[n], y[n])
          arrow_head = (numpy.mean(x[n:n + 2]), numpy.mean(y[n:n + 2]))
          p = mpatches.FancyArrowPatch(
            arrow_tail, arrow_head, transform=transform,
            **arrow_kw)
          axes.add_patch(p)
          arrows.append(p)
     return arrows