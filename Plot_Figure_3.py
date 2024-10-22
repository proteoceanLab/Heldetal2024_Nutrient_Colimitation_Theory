import colimitation_models
import colimitation_plots
import matplotlib
from matplotlib import pyplot
import numpy
import pandas
from matplotlib.ticker import LogLocator


def main():

     # Fraction of limitation in a single resource to define single limitation
     single_lim_threshold = 0.95

     # Fraction of limitation in a single resource to define double colimitation
     #double_colim_threshold = 0.5

     # Set up log meshes for Meff contours
     R1_list = numpy.logspace(-4, 4, 300)
     R2_list = numpy.logspace(-4.5, 4.5, 300)
     R1_mesh, R2_mesh = numpy.meshgrid(R1_list, R2_list)

     # Read model parameters for growth rate
     rate_model_to_plot = "pat_rmin"
     dataframe = pandas.read_csv("./Data/Rate_data_fits_rep123.csv", index_col=False)
     dataframe_filtered = dataframe[dataframe["Model name"] == rate_model_to_plot]
     row_index = dataframe_filtered.index[0]
     gmax, a1, a2, R1min, R2min = [float(x) for x in str(dataframe_filtered["Parameter values"][row_index]).split(";")]

     # Iterate over points in R1-R2 space and calculate limitation coefficients and Meff
     L1_rate_mesh = numpy.zeros(R1_mesh.shape)
     L2_rate_mesh = numpy.zeros(R1_mesh.shape)
     Meff_rate_mesh = numpy.zeros(R1_mesh.shape)
     for i in range(len(R2_list)):
          for j in range(len(R1_list)):
               L1_rate_mesh[i][j] = colimitation_models.CalcLimCoeff(0, [R1_list[j], R2_list[i]], [a1, a2], gmax, rate_model_to_plot.strip("_rmin"))
               L2_rate_mesh[i][j] = colimitation_models.CalcLimCoeff(1, [R1_list[j], R2_list[i]], [a1, a2], gmax, rate_model_to_plot.strip("_rmin"))
               Meff_rate_mesh[i][j] = colimitation_models.CalcMeff([R1_list[j], R2_list[i]], [a1, sa], gmax, rate_model_to_plot.strip("_rmin"))

     # Read model parameters for growth yield
     yield_model_to_plot = "pat"
     dataframe = pandas.read_csv("./Data/Yield_data_fits_rep123.csv", index_col=False)
     dataframe_filtered = dataframe[dataframe["Model name"] == yield_model_to_plot]
     row_index = dataframe_filtered.index[0]
     Nmax, s1, s2 = [float(x) for x in str(dataframe_filtered["Parameter values"][row_index]).split(";")]

     # Iterate over points in R1-R2 space and calculate limitation coefficients and Meff
     L1_yield_mesh = numpy.zeros(R1_mesh.shape)
     L2_yield_mesh = numpy.zeros(R1_mesh.shape)
     Meff_yield_mesh = numpy.zeros(R1_mesh.shape)
     for i in range(len(R2_list)):
          for j in range(len(R1_list)):
               L1_yield_mesh[i][j] = colimitation_models.CalcLimCoeff(0, [R1_list[j], R2_list[i]], [s1, s2], Nmax, yield_model_to_plot.strip("_rmin"))
               L2_yield_mesh[i][j] = colimitation_models.CalcLimCoeff(1, [R1_list[j], R2_list[i]], [s1, s2], Nmax, yield_model_to_plot.strip("_rmin"))
               Meff_yield_mesh[i][j] = colimitation_models.CalcMeff([R1_list[j], R2_list[i]], [s1, s2], Nmax, yield_model_to_plot.strip("_rmin"))

################################################################################

     # Get rate and yield model parameters for all replicates and bootstraps

     # Number of experimental replicates and number of bootstrapped data sets
     # from those replicates
     num_replicates = 3
     num_bootstraps = 100

     # Colors for the experimental replicates
     rep_colors = ["tab:blue", "tab:red", "tab:orange"]

     # Lists of fitted parameters (glucose-ammonium stoichiometry, glucose threshold,
     # ammonium threshold)
     rate_params_all_reps = []
     rate_params_indiv_reps = []
     rate_params_bootstraps = []

     # Get dataframes for fit to all replicates, to each replicate individually,
     # and to all bootstrapped data sets.  Shape is (num replicates or bootstraps,
     # num models, num parameters).  
     rate_data_all_reps = pandas.read_csv("./Data/Rate_data_fits_rep123.csv", index_col=False)
     rate_data_indiv_reps = [pandas.read_csv("./Data/Rate_data_fits_rep" + str(r + 1) + ".csv", index_col=False) for r in range(num_replicates)]
     rate_data_bootstraps = pandas.read_csv("./Data/Rate_bootstraps_fits.csv", index_col=False)

     # For this model get parameters from fits to all replicates
     df_this_model = rate_data_all_reps[rate_data_all_reps["Model name"] == rate_model_to_plot]
     row_index = df_this_model.index[0]
     params = [float(x) for x in str(df_this_model["Parameter values"][row_index]).split(";")]
     gmax, s1, s2 = params[:3]
     rate_params_all_reps = [s2/s1, gmax/s1, gmax/s2]

     # For this model get parameters from fits to individual replicates
     for r in range(num_replicates):
          df_this_model = rate_data_indiv_reps[r][rate_data_indiv_reps[r]["Model name"] == rate_model_to_plot]
          row_index = df_this_model.index[0]
          params = [float(x) for x in str(df_this_model["Parameter values"][row_index]).split(";")]
          gmax, a1, a2 = params[:3]
          rate_params_indiv_reps.append([a2/a1, gmax/a1, gmax/a2])

     # For this model get parameters from fits to each bootstrapped data set
     for b in range(num_bootstraps):
          df_this_model = rate_data_bootstraps[rate_data_bootstraps["Model name"] == rate_model_to_plot]
          row_index = df_this_model.index[0]
          params = []
          for x in str(df_this_model[f"Parameter values bootstrap {b + 1}"][row_index]).split(";"):
               if x == "NA":
                    params.append(numpy.nan)
               else:
                    params.append(float(x))
          gmax, a1, a2 = params[:3]
          rate_params_bootstraps.append([a2/a1, gmax/a1, gmax/a2])

     # Lists for fitted parameters
     yield_params_all_reps = []
     yield_params_indiv_reps = []
     yield_params_bootstraps = []

     # Get dataframes for fit to all replicates, to each replicate individually,
     # and to all bootstrapped data sets.  Shape is (num replicates or bootstraps,
     # num models, num parameters).  
     yield_data_all_reps = pandas.read_csv("./Data/Yield_data_fits_rep123.csv", index_col=False)
     yield_data_indiv_reps = [pandas.read_csv("./Data/Yield_data_fits_rep" + str(r + 1) + ".csv", index_col=False) for r in range(num_replicates)]
     yield_data_bootstraps = pandas.read_csv("./Data/Yield_bootstraps_fits.csv", index_col=False)

     # For this model get parameters from fits to all replicates
     df_this_model = yield_data_all_reps[yield_data_all_reps["Model name"] == yield_model_to_plot]
     row_index = df_this_model.index[0]
     params = [float(x) for x in str(df_this_model["Parameter values"][row_index]).split(";")]
     Nmax, s1, s2 = params[:3]
     yield_params_all_reps = [s2/s1, Nmax/s1, Nmax/s2]

     # For this model get parameters from fits to individual replicates
     for r in range(num_replicates):
          df_this_model = yield_data_indiv_reps[r][yield_data_indiv_reps[r]["Model name"] == yield_model_to_plot]
          row_index = df_this_model.index[0]
          params = [float(x) for x in str(df_this_model["Parameter values"][row_index]).split(";")]
          Nmax, s1, s2 = params[:3]
          yield_params_indiv_reps.append([s2/s1, Nmax/s1, Nmax/s2])

     # For this model get parameters from fits to each bootstrapped data set
     for b in range(num_bootstraps):
          df_this_model = yield_data_bootstraps[yield_data_bootstraps["Model name"] == yield_model_to_plot]
          row_index = df_this_model.index[0]
          params = []
          for x in str(df_this_model[f"Parameter values bootstrap {b + 1}"][row_index]).split(";"):
               if x == "NA":
                    params.append(numpy.nan)
               else:
                    params.append(float(x))
          Nmax, s1, s2 = params[:3]
          yield_params_bootstraps.append([s2/s1, Nmax/s1, Nmax/s2])

################################################################################

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Initialize figure
     figure = pyplot.figure(figsize=(3*4, 3*1))
     figure.subplots_adjust(wspace=-0.25, hspace=0.1)

     # Set column widths, pads, and coordinates
     column_widths = [1, 2, 1, 1]
     column_pads = [1, 1, 1, 0]
     column_coords = [sum(column_widths[:i]) + sum(column_pads[:i]) for i in range(len(column_widths))]

################################################################################

     # Plot contours of rate limitation coefficients
     axis = pyplot.subplot2grid((2, column_coords[-1] + 1), (0, column_coords[0]), colspan=column_widths[0])
     axis.set_ylabel("Ammonium (mM)", fontsize=colimitation_plots.axis_label_size)
     axis.set_xscale("log")
     axis.set_yscale("log")
     axis.set_xticks(numpy.logspace(-4, 4, 5))
     axis.set_yticks(numpy.logspace(-4, 4, 5))
     axis.set_ylim([min(R2_list), max(R2_list)])
     axis.set_aspect("equal")
     pyplot.setp(axis.get_xticklabels(), visible=False)

     limitation_levels = [1e-4, 0.01, 0.9, 1]
     cmap = matplotlib.colormaps["Blues"]
     limitation_color_levels_1 = [cmap(c) for c in numpy.linspace(0.4, 1, len(limitation_levels))]
     cmap = matplotlib.colormaps["Oranges"]
     limitation_color_levels_2 = [cmap(c) for c in numpy.linspace(0.4, 1, len(limitation_levels))]
     cmap = matplotlib.colormaps["Greens"]
     limitation_color_levels_3 = [cmap(c) for c in numpy.linspace(0.4, 1, len(limitation_levels))]
     contours = axis.contourf(R1_mesh, R2_mesh, L1_rate_mesh, alpha=0.5, colors=limitation_color_levels_1, levels=limitation_levels) 
     contours = axis.contourf(R1_mesh, R2_mesh, L2_rate_mesh, alpha=0.5, colors=limitation_color_levels_2, levels=limitation_levels)

     axis.annotate(r"$L^\mathrm{rate}_\mathrm{glu}$", xytext=(1e1, 10**(3.5)), xy=(1e-4, 10**(3.5)), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.tick_label_size, color="black"), color="black", ha="center", va="center", size=colimitation_plots.legend_label_size)
     axis.annotate(r"$L^\mathrm{rate}_\mathrm{amm}$", xytext=(10**(2.75), 1e1), xy=(10**(2.75), 1e-4), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.tick_label_size, color="black"), color="black", ha="center", va="center", size=colimitation_plots.legend_label_size)

################################################################################

     # Plot contours of yield limitation coefficients
     axis = pyplot.subplot2grid((2, column_coords[-1] + 1), (1, column_coords[0]), colspan=column_widths[0])
     axis.set_xlabel("Glucose (mM)", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("Ammonium (mM)", fontsize=colimitation_plots.axis_label_size)
     axis.set_xscale("log")
     axis.set_yscale("log")
     axis.set_xticks(numpy.logspace(-4, 4, 5))
     axis.tick_params(axis='x', which='major', labelrotation=45)
     axis.set_yticks(numpy.logspace(-4, 4, 5))
     axis.set_ylim([min(R2_list), max(R2_list)])
     axis.set_aspect("equal")

     limitation_levels = [0.01, 0.9, 0.999, 1]
     cmap = matplotlib.colormaps["Blues"]
     limitation_color_levels_1 = [cmap(c) for c in numpy.linspace(0.4, 1, len(limitation_levels))]
     cmap = matplotlib.colormaps["Oranges"]
     limitation_color_levels_2 = [cmap(c) for c in numpy.linspace(0.4, 1, len(limitation_levels))]
     cmap = matplotlib.colormaps["Greens"]
     limitation_color_levels_3 = [cmap(c) for c in numpy.linspace(0.4, 1, len(limitation_levels))]
     contours = axis.contourf(R1_mesh, R2_mesh, L1_yield_mesh, alpha=0.5, colors=limitation_color_levels_1, levels=limitation_levels) 
     contours = axis.contourf(R1_mesh, R2_mesh, L2_yield_mesh, alpha=0.5, colors=limitation_color_levels_2, levels=limitation_levels)

     axis.annotate(r"$L^\mathrm{yield}_\mathrm{glu}$", xytext=(1e1, 10**(3.5)), xy=(1e-4, 10**(3.5)), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.tick_label_size, color="black"), color="black", ha="center", va="center", size=colimitation_plots.legend_label_size)
     axis.annotate(r"$L^\mathrm{yield}_\mathrm{amm}$", xytext=(10**(2.75), 1e1), xy=(10**(2.75), 1e-4), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.tick_label_size, color="black"), color="black", ha="center", va="center", size=colimitation_plots.legend_label_size)

################################################################################

     # Phase diagram of rate and yield colimitation regimes
     axis = pyplot.subplot2grid((2, column_coords[-1] + 1), (0, column_coords[1]), colspan=column_widths[1], rowspan=2)
     axis.text(-1.1, 1.05, "A", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.text(-0.25, 1.05, "B", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_xlabel("Glucose concentration (mM)", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("Ammonium concentration (mM)", fontsize=colimitation_plots.axis_label_size)
     axis.set_aspect("equal")

     # Draw shaded regions for rate and yield colimitation
     contourfs = axis.contourf(R1_mesh, R2_mesh, Meff_rate_mesh, colors=["blue"], levels=[1/single_lim_threshold, 3], zorder=0, alpha=0.2)
     #contourfs = axis.contourf(R1_mesh, R2_mesh, Meff_rate_mesh, colors=["blue"], levels=[1/double_colim_threshold, 3], zorder=0, alpha=1)
     contourfs = axis.contourf(R1_mesh, R2_mesh, Meff_yield_mesh, colors=["red"], levels=[1/single_lim_threshold, 3], zorder=1, alpha=0.2)
     #contourfs = axis.contourf(R1_mesh, R2_mesh, Meff_yield_mesh, colors=["red"], levels=[1/double_colim_threshold, 3], zorder=1, alpha=1)
     axis.text(8000, 2, r"\textbf{Yield}" + "\n" + r"\textbf{colim}", color="red", ha="right", va="center", fontsize=colimitation_plots.axis_label_size)
     axis.text(3e-2, 7000, r"\textbf{Rate}" + "\n" + r"\textbf{colim}", color="blue", ha="center", va="top", fontsize=colimitation_plots.axis_label_size)

     # Add box for resource concentration range in natural environments
     glucose_min, glucose_max = 0.4, 24
     ammonium_min, ammonium_max = 0.2, 0.31
     axis.fill_between([glucose_min, glucose_max], [ammonium_min, ammonium_min], y2=ammonium_max, color="0.4", alpha=0.75, zorder=2)
     axis.annotate(r"\textbf{Natural}" + "\n" + r"\textbf{environments}", xytext=(9e-1, 3e-3), xy=(numpy.sqrt(glucose_min*glucose_max), ammonium_min), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.tick_label_size, color="0.4"), color="0.4", ha="center", va="center", size=colimitation_plots.legend_label_size)

     # Typical laboratory concentrations
     glucose_point2percent = 11.101489819933835
     ammonium_M9 = 18.7
     axis.scatter(glucose_point2percent, ammonium_M9, marker="o", color="black", zorder=3)
     axis.annotate(r"\textbf{Typical M9}" + "\n" + r"\textbf{0.2\% glucose}", xytext=(3e2, 1e3), xy=(glucose_point2percent, ammonium_M9), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.tick_label_size, color="black"), ha="center", va="center", size=colimitation_plots.legend_label_size)

     # LTEE concentrations (but in Davis-Mingioli not M9)
     #glucose_LTEE = 0.13876862274917295
     #ammonium_LTEE = 15.2
     #axis.scatter(glucose_LTEE, ammonium_LTEE, marker="s", color="tab:green", zorder=3) 
     #axis.annotate(r"\textbf{LTEE}", xytext=(2e-2, 2), xy=(glucose_LTEE, ammonium_LTEE), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.tick_label_size, color="tab:green"), color="tab:green", ha="center", va="center", size=colimitation_plots.legend_label_size)

     # Concentrations used in Bren et al. 2013 BMC Syst Biol as single limitation 
     # for glucose or ammonium (combined with standard lab concentrations)
     glucose_Bren = 0.14
     ammonium_Bren = 0.24
     axis.scatter(glucose_Bren, ammonium_M9, marker="^", s=20, color="tab:green", zorder=3)
     axis.scatter(glucose_point2percent, ammonium_Bren, marker="+", s=50, color="tab:pink", zorder=3)
     axis.annotate(r"\textbf{LTEE/}" + "\n" + r"\textbf{Bren et al.}" + "\n" + r"\textbf{Glu single}" + "\n" + r"\textbf{lim}", xytext=(2e-3, 5e0), xy=(glucose_Bren, ammonium_M9), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.tick_label_size, color="tab:green"), color="tab:green", ha="center", va="top", size=colimitation_plots.legend_label_size)
     axis.annotate(r"\textbf{Bren}" + "\n" + r"\textbf{et al.}" + "\n" + r"\textbf{Amm}" + "\n" + r"\textbf{single}" + "\n" + r"\textbf{lim}", xytext=(2e3, 1e-1), xy=(glucose_point2percent, ammonium_Bren), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.tick_label_size, color="tab:pink"), color="tab:pink", ha="center", va="top", size=colimitation_plots.legend_label_size)

     # Concentration at which ammonium transporter turns on (from Kim et al. 
     # 2012 Mol Syst Biol)
     glucose_Kim = 2*glucose_point2percent
     ammonium_Kim = 0.03
     axis.scatter(glucose_Kim, ammonium_Kim, marker="*", color="tab:brown", zorder=3)   
     axis.annotate(r"\textbf{Kim et al.}" + "\n" + r"\textbf{Amm transport}", xytext=(1e2, 1e-3), xy=(glucose_Kim, ammonium_Kim), arrowprops=dict(arrowstyle="-|>", mutation_scale=colimitation_plots.tick_label_size, color="tab:brown"), color="tab:brown", ha="center", va="top", size=colimitation_plots.legend_label_size)

     # Set axis ranges, scales, and ticks (force log minor ticks for x-axis)
     axis.set_xscale("log")
     axis.set_yscale("log")
     axis.set_xlim([1e-4, 1e4])
     axis.set_ylim([1e-4, 1e4])
     axis.set_xticks(numpy.logspace(-4, 4, 9))
     axis.set_yticks(numpy.logspace(-4, 4, 9))
     axis.tick_params(axis="x", which="major", labelrotation=45)
     axis.xaxis.get_minor_locator().set_params(numticks=99, subs=numpy.linspace(0.1, 0.9, 9))

################################################################################

     # Plot glucose-ammonium stoichiometry
     axis = pyplot.subplot2grid((2, column_coords[-1] + 1), (0, column_coords[2]), colspan=column_widths[2], rowspan=2)
     axis.text(-0.32, 1.05, "C", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.text(1.1, 1.05, "D", transform=axis.transAxes, fontsize=colimitation_plots.panel_label_size)
     axis.set_ylabel("Glucose-ammonium stoichiometry\nunder colimitation", fontsize=colimitation_plots.axis_label_size)
     axis.set_xlim([-0.5, 1.5])
     axis.set_ylim([1, 7])

     # Plot points and bootstrap boxes for rate parameters
     axis.scatter(0, rate_params_all_reps[0], color="black", marker="o", s=40, zorder=1)
     for r in range(num_replicates):
          axis.scatter(0, rate_params_indiv_reps[r][0], color=rep_colors[r], marker="_", s=300, zorder=2)
     params_bootstraps = numpy.array([rate_params_bootstraps[b][0] for b in range(num_bootstraps)])
     boxes = axis.boxplot(params_bootstraps, positions=[0], sym="", patch_artist=True, zorder=0, 
          medianprops={"color": "0.8", "linewidth": 0.5},
          boxprops={"facecolor": "0.8", "edgecolor": "0.8", "linewidth": 10},
          whiskerprops={"color": "0.8", "linewidth": 1.5},
          capprops={"color": "0.8", "linewidth": 1.5})

     # Plot points and bootstrap boxes for yield parameters
     axis.scatter(1, yield_params_all_reps[0], color="black", marker="o", s=40, zorder=1)
     for r in range(num_replicates):
          axis.scatter(1, yield_params_indiv_reps[r][0], color=rep_colors[r], marker="_", s=300, zorder=2)
     params_bootstraps = numpy.array([yield_params_bootstraps[b][0] for b in range(num_bootstraps)])
     boxes = axis.boxplot(params_bootstraps, positions=[1], sym="", patch_artist=True, zorder=0, 
          medianprops={"color": "0.8", "linewidth": 0.5},
          boxprops={"facecolor": "0.8", "edgecolor": "0.8", "linewidth": 10},
          whiskerprops={"color": "0.8", "linewidth": 1.5},
          capprops={"color": "0.8", "linewidth": 1.5})

     # Plot dummy line with the same color as boxes but so it appears in 
     # reasonable size in the legend
     axis.scatter([], [], color="black", marker="o", s=20, label="Combined\nreps")
     for r in range(num_replicates):
          axis.plot([], [], "-", color=rep_colors[r], linewidth=2, label=f"Rep {r + 1}")
     axis.plot([], [], "-", color="0.8", linewidth=5, label="Boot-\nstrapped\ndata")

     axis.set_xticks([0, 1], labels=["Growth\nrate", "Growth\nyield"], fontsize=colimitation_plots.axis_label_size)
     axis.legend(loc="upper right", fontsize=0.9*colimitation_plots.legend_label_size, handlelength=0.5)

################################################################################

     # Plot glucose concentration at colimitation with implicit factor
     axis = pyplot.subplot2grid((2, column_coords[-1] + 1), (0, column_coords[3]), colspan=column_widths[3])
     axis.set_ylabel("Glucose (mM)", fontsize=colimitation_plots.axis_label_size)
     axis.set_xlim([-0.5, 1.5])
     axis.set_ylim([10**(-2.25), 1e1])
     axis.set_yscale("log")
     axis.set_title("Concentration at colimitation\nwith implicit factors", fontsize=colimitation_plots.axis_label_size)

     # Plot points and bootstrap boxes for rate parameters
     axis.scatter(0, rate_params_all_reps[1], color="black", marker="o", s=40, zorder=1)
     for r in range(num_replicates):
          axis.scatter(0, rate_params_indiv_reps[r][1], color=rep_colors[r], marker="_", s=300, zorder=2)
     params_bootstraps = numpy.array([rate_params_bootstraps[b][1] for b in range(num_bootstraps)])
     axis.boxplot(params_bootstraps, positions=[0], sym="", patch_artist=True, zorder=0, 
          medianprops={"color": "0.8", "linewidth": 0.5},
          boxprops={"facecolor": "0.8", "edgecolor": "0.8", "linewidth": 10},
          whiskerprops={"color": "0.8", "linewidth": 1.5},
          capprops={"color": "0.8", "linewidth": 1.5})

     # Plot points and bootstrap boxes for yield parameters
     axis.scatter(1, yield_params_all_reps[1], color="black", marker="o", s=40, zorder=1)
     for r in range(num_replicates):
          axis.scatter(1, yield_params_indiv_reps[r][1], color=rep_colors[r], marker="_", s=300, zorder=2)
     params_bootstraps = numpy.array([yield_params_bootstraps[b][1] for b in range(num_bootstraps)])
     axis.boxplot(params_bootstraps, positions=[1], sym="", patch_artist=True, zorder=0, 
          medianprops={"color": "0.8", "linewidth": 0.5},
          boxprops={"facecolor": "0.8", "edgecolor": "0.8", "linewidth": 10},
          whiskerprops={"color": "0.8", "linewidth": 1.5},
          capprops={"color": "0.8", "linewidth": 1.5})

     # Set y-ticks and force log minor ticks
     axis.set_yticks(numpy.logspace(-2, 1, 4))
     axis.yaxis.get_minor_locator().set_params(numticks=99, subs=numpy.linspace(0.1, 0.9, 9))

################################################################################

     # Plot ammonium concentration at colimitation with implicit factor
     axis = pyplot.subplot2grid((2, column_coords[-1] + 1), (1, column_coords[3]), colspan=column_widths[3])
     axis.set_ylabel("Ammonium (mM)", fontsize=colimitation_plots.axis_label_size)
     axis.set_xlim([-0.5, 1.5])
     axis.set_ylim([1e-3, 10**(1.25)])
     axis.set_yscale("log")

     # Plot points and bootstrap boxes for rate parameters
     axis.scatter(0, rate_params_all_reps[2], color="black", marker="o", s=40, zorder=1)
     for r in range(num_replicates):
          axis.scatter(0, rate_params_indiv_reps[r][2], color=rep_colors[r], marker="_", s=300, zorder=2)
     params_bootstraps = numpy.array([rate_params_bootstraps[b][2] for b in range(num_bootstraps)])
     axis.boxplot(params_bootstraps, positions=[0], sym="", patch_artist=True, zorder=0, 
          medianprops={"color": "0.8", "linewidth": 0.5},
          boxprops={"facecolor": "0.8", "edgecolor": "0.8", "linewidth": 10},
          whiskerprops={"color": "0.8", "linewidth": 1.5},
          capprops={"color": "0.8", "linewidth": 1.5})

     # Plot points and bootstrap boxes for yield parameters
     axis.scatter(1, yield_params_all_reps[2], color="black", marker="o", s=40, zorder=1)
     for r in range(num_replicates):
          axis.scatter(1, yield_params_indiv_reps[r][2], color=rep_colors[r], marker="_", s=300, zorder=2)
     params_bootstraps = numpy.array([yield_params_bootstraps[b][2] for b in range(num_bootstraps)])
     axis.boxplot(params_bootstraps, positions=[1], sym="", patch_artist=True, zorder=0, 
          medianprops={"color": "0.8", "linewidth": 0.5},
          boxprops={"facecolor": "0.8", "edgecolor": "0.8", "linewidth": 10},
          whiskerprops={"color": "0.8", "linewidth": 1.5},
          capprops={"color": "0.8", "linewidth": 1.5})

     # Set y-ticks and force log minor ticks
     axis.set_xticks([0, 1], labels=["Growth\nrate", "Growth\nyield"], fontsize=colimitation_plots.axis_label_size)
     axis.set_yticks(numpy.logspace(-3, 1, 5))
     axis.yaxis.get_minor_locator().set_params(numticks=99, subs=numpy.linspace(0.1, 0.9, 9))

################################################################################

     figure.savefig("Figure_3.pdf", bbox_inches="tight")

################################################################################


if __name__ == '__main__':
     main()
