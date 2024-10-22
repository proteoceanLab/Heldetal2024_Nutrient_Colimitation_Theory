import colimitation_data_analysis
import colimitation_models
import colimitation_plots
from matplotlib import pyplot
import pandas
import numpy
from scipy import stats


def main():    

     # Number of experimental replicates and number of bootstrapped data sets
     # from those replicates
     num_replicates = 3
     num_bootstraps = 100

     # Colors for the experimental replicates
     rep_colors = ["tab:blue", "tab:red", "tab:orange"]

     # Lists for fitted parameters
     fit_params_all_reps = []
     fit_params_indiv_reps = [[] for r in range(num_replicates)]
     fit_params_bootstraps = [[] for b in range(num_bootstraps)]

     # Get dataframes for fit to all replicates, to each replicate individually,
     # and to all bootstrapped data sets.  Shape is (num replicates or bootstraps,
     # num models, num parameters).  
     fit_data_all_reps = pandas.read_csv("./Data/Rate_data_fits_rep123.csv", index_col=False)
     fit_data_indiv_reps = [pandas.read_csv("./Data/Rate_data_fits_rep" + str(r + 1) + ".csv", index_col=False) for r in range(num_replicates)]
     fit_data_bootstraps = pandas.read_csv("./Data/Rate_bootstraps_fits.csv", index_col=False)

     # Get list of indices for models that have Rmin
     rmin_model_indices = [i for i in range(len(colimitation_models.all_2D_trait_models)) if "rmin" in colimitation_models.all_2D_trait_models[i]]

     # Iterate over all fitted models
     for m in range(len(colimitation_models.all_2D_trait_models)):

          model_name = colimitation_models.all_2D_trait_models[m]

          # For this model get parameters from fits to all replicates
          df_this_model = fit_data_all_reps[fit_data_all_reps["Model name"] == model_name]
          row_index = df_this_model.index[0]
          Rsq = float(df_this_model["R^2"])
          akaike_weight = float(df_this_model["Akaike weight"])
          params = [float(x) for x in str(df_this_model["Parameter values"][row_index]).split(";")]
          zmax, s1, s2 = params[:3]
          if "rmin" not in model_name:
               fit_data = [Rsq, akaike_weight, zmax, s1, s2, s2/s1, zmax/s1, zmax/s2]
          else:      
               R1min, R2min = params[-2:]
               fit_data = [Rsq, akaike_weight, zmax, s1, s2, s2/s1, zmax/s1, zmax/s2, R1min, R2min]
          fit_params_all_reps.append(fit_data)     

          # For this model get parameters from fits to individual replicates
          for r in range(num_replicates):
               df_this_model = fit_data_indiv_reps[r][fit_data_indiv_reps[r]["Model name"] == model_name]
               row_index = df_this_model.index[0]
               Rsq = float(df_this_model["R^2"])
               akaike_weight = float(df_this_model["Akaike weight"])
               params = [float(x) for x in str(df_this_model["Parameter values"][row_index]).split(";")]
               zmax, s1, s2 = params[:3]
               if "rmin" not in model_name:
                    fit_data = [Rsq, akaike_weight, zmax, s1, s2, s2/s1, zmax/s1, zmax/s2]
               else:      
                    R1min, R2min = params[-2:]
                    fit_data = [Rsq, akaike_weight, zmax, s1, s2, s2/s1, zmax/s1, zmax/s2, R1min, R2min]
               fit_params_indiv_reps[r].append(fit_data)

          # For this model get parameters from fits to each bootstrapped data set
          for b in range(num_bootstraps):
               df_this_model = fit_data_bootstraps[fit_data_bootstraps["Model name"] == model_name]
               row_index = df_this_model.index[0]
               Rsq = float(df_this_model[f"R^2 bootstrap {b + 1}"])
               akaike_weight = float(df_this_model[f"Akaike weight bootstrap {b + 1}"])
               params = []
               for x in str(df_this_model[f"Parameter values bootstrap {b + 1}"][row_index]).split(";"):
                    if x == "NA":
                         params.append(numpy.nan)
                    else:
                         params.append(float(x))
               zmax, s1, s2 = params[:3]
               if "rmin" not in model_name:
                    fit_data = [Rsq, akaike_weight, zmax, s1, s2, s2/s1, zmax/s1, zmax/s2]
               else:      
                    R1min, R2min = params[-2:]
                    fit_data = [Rsq, akaike_weight, zmax, s1, s2, s2/s1, zmax/s1, zmax/s2, R1min, R2min]
               fit_params_bootstraps[b].append(fit_data)

     # Axis labels for each fit parameter to plot
     ylabels = [    "Quality of growth rate\nmodel fit $R^2$", 
                    "Relative Akaike weight", 
                    "Maximum growth rate\n$g_\mathrm{max}$ (per hour)",
                    "Glucose affinity $a_\mathrm{glu}$\n(per hour per mM glucose)",
                    "Ammonium affinity $a_\mathrm{amm}$\n(per hour per mM ammonium)",
                    "Glucose-ammonium\nstoichiometry $a_\mathrm{amm}/a_\mathrm{glu}$\n(mM glucose per mM ammonium)",
                    "Glucose half-saturation\nconcentration $g_\mathrm{max}/a_\mathrm{glu}$ (mM)",
                    "Ammonium half-saturation\nconcentration $g_\mathrm{max}/a_\mathrm{amm}$ (mM)",
                    "Minimum glucose\nconcentration $R_\mathrm{glu,min}$ (mM)",
                    "Minimum ammonium\nconcentration $R_\mathrm{amm,min}$ (mM)"]

     # Axis limits for each fit parameter to plot
     ylims = [ (0, 1),
               (10**(-12), 10**(0)),
               (0.5, 0.9),
               (0, 200),
               (0, 200),
               (0, 8),
               (0, 0.05),
               (0, 0.015),
               (0, 0.1),
               (0, 0.05)]

     # Perform Mann-Whitney U test to compare R^2 distributions from the bootstrapped
     # data sets, just between the top few models in the all-replicate fits
     print("p-values for Mann-Whitney U test between bootstrapped R^2 distributions")
     ref_model_index = 0
     for m in range(len(colimitation_models.all_2D_trait_models)):
          if m == ref_model_index: continue
          Rsq1 = [fit_params_bootstraps[b][ref_model_index][0] for b in range(num_bootstraps)]
          Rsq2 = [fit_params_bootstraps[b][m][0] for b in range(num_bootstraps)]
          result = stats.mannwhitneyu(Rsq1, Rsq2)
          print("\t", colimitation_models.all_2D_trait_models[ref_model_index], "vs.", colimitation_models.all_2D_trait_models[m], "p =", result.pvalue)

################################################################################

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Initialize figure
     figure = pyplot.figure(figsize=(18, 3*len(ylabels)/2))
     figure.subplots_adjust(hspace=0.1, wspace=0.15)

     # Iterate over each fit parameter to plot
     for p in range(len(ylabels)):

          # Set up axis
          axis = figure.add_subplot(int(len(ylabels)/2), 2, p + 1)
          axis.set_ylabel(ylabels[p], fontsize=colimitation_plots.axis_label_size)
          axis.set_xlim([-0.5, len(colimitation_models.all_2D_trait_models) - 0.5])
          axis.set_ylim(ylims[p])
          if "Akaike" in ylabels[p]:
               axis.set_yscale("log")
          for i in range(len(colimitation_models.all_2D_trait_models)):
               axis.axvline(i, linestyle="--", linewidth=0.25, color="0.8", zorder=-1)

          # If the parameter to plot is Rmin, use only a subset of model indices
          # for Rmin models
          if "Minimum" in ylabels[p]:
               model_indices = rmin_model_indices
          else:
               model_indices = range(len(colimitation_models.all_2D_trait_models))
     
          # Plot fit parameters for all-replicate fit
          axis.scatter(model_indices, [fit_params_all_reps[i][p] for i in model_indices], color="black", s=10, zorder=1, label="Combined reps")
     
          # Plot fit parameters for individual replicates
          for r in range(num_replicates):
               axis.scatter(model_indices, [fit_params_indiv_reps[r][i][p] for i in model_indices], marker="_", label="Rep " + str(r + 1), color=rep_colors[r], zorder=2)
     
          # Plot fit parameters for bootstrapped data
          params_bootstraps = numpy.array([[fit_params_bootstraps[b][i][p] for i in model_indices] for b in range(num_bootstraps)])
          boxes = axis.boxplot(params_bootstraps, positions=model_indices, sym="", patch_artist=True, zorder=0, 
               medianprops={"color": "0.8", "linewidth": 0.5},
               boxprops={"facecolor": "0.8", "edgecolor": "0.8", "linewidth": 0.5},
               whiskerprops={"color": "0.8", "linewidth": 1.5},
               capprops={"color": "0.8", "linewidth": 1.5}, label="Bootstrapped data")

          # Show legend in first plot only
          if p == 0:
               axis.legend(loc=(0, 1.8), fontsize=colimitation_plots.axis_label_size)

          # Add top tick marks for first two panels
          if p == 0 or p == 1:
               axis_twiny = axis.twiny()
               axis_twiny.set_xticks(range(len(colimitation_models.all_2D_trait_models)), labels=colimitation_models.all_2D_trait_models_formatted, fontsize=colimitation_plots.tick_label_size, rotation=90)
               axis_twiny.set_xlim([-0.5, len(colimitation_models.all_2D_trait_models) - 0.5])

          # Label bottom ticks in last two plots only
          if p == (len(ylabels) - 2) or p == (len(ylabels) - 1):
               axis.set_xticks(range(len(colimitation_models.all_2D_trait_models)), labels=colimitation_models.all_2D_trait_models_formatted, fontsize=colimitation_plots.tick_label_size, rotation=90)
          else:
               axis.set_xticks(range(len(colimitation_models.all_2D_trait_models)), labels=[])

     figure.savefig("Figure_S14.pdf", bbox_inches="tight")

################################################################################


if __name__ == '__main__':
     main()
