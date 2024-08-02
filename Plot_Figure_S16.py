import colimitation_data_analysis
import colimitation_models
import colimitation_plots
from matplotlib import pyplot
import pandas
import numpy
from scipy import stats


def main():    

     # Get true parameter values that were used for simulations
     dataframe = pandas.read_csv("./Data/Rate_data_fits_rep123.csv", index_col=False)
     true_model = "liebig_blackman_rmin"
     true_model_formatted = colimitation_models.all_2D_trait_models_formatted[colimitation_models.all_2D_trait_models.index(true_model)]
     dataframe_filtered = dataframe[dataframe["Model name"] == true_model]
     row_index = dataframe_filtered.index[0]
     zmax, s1, s2, R1min, R2min = [float(x) for x in str(dataframe_filtered["Parameter values"][row_index]).split(";")]
     true_parameters = [zmax, s1, s2, s2/s1, zmax/s1, zmax/s2, R1min, R2min]

     # Number of simulated data sets
     num_sims = 10000

     # Lists for fitted parameters
     fit_params_sims = [[] for s in range(num_sims)]

     # Get list of indices for models that have Rmin
     rmin_model_indices = [i for i in range(len(colimitation_models.all_2D_trait_models)) if "rmin" in colimitation_models.all_2D_trait_models[i]]

     # Read fit data for simulations
     fit_data_sim = pandas.read_csv(f"./Data/Rate_sim_{true_model}_fits.csv", index_col=False)

     # Get parameters from fits to each simulated data set
     for s in range(num_sims):

          # Iterate over all fitted models
          for model_name in colimitation_models.all_2D_trait_models:

               df_this_model = fit_data_sim[fit_data_sim["Model name"] == model_name]
               row_index = df_this_model.index[0]
               Rsq = float(df_this_model[f"R^2 sim {s + 1}"])
               akaike_weight = float(df_this_model[f"Akaike weight sim {s + 1}"])
               params = []
               for x in str(df_this_model[f"Parameter values sim {s + 1}"][row_index]).split(";"):
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
               fit_params_sims[s].append(fit_data)

     # Axis labels for each fit parameter to plot
     ylabels = [    "Quality of growth rate\nmodel fit $R^2$", 
                    "Relative Akaike weight", 
                    "Maximum growth rate\n$g_\mathrm{max}$ (per hour)",
                    "Growth rate per\nglucose $a_\mathrm{glucose}$\n(per hour per mM glucose)",
                    "Growth rate per\nammonium $a_\mathrm{ammonium}$\n(per hour per mM ammonium)",
                    "Glucose-ammonium\nstoichiometry $a_\mathrm{ammonium}/a_\mathrm{glucose}$\n(mM glucose per mM ammonium)",
                    "Glucose threshold\n$g_\mathrm{max}/a_\mathrm{glucose}$ (mM)",
                    "Ammonium threshold\n$g_\mathrm{max}/a_\mathrm{ammonium}$ (mM)",
                    "Minimum glucose\nconcentration $R_\mathrm{glucose,min}$ (mM)",
                    "Minimum ammonium\nconcentration $R_\mathrm{ammonium,min}$ (mM)"]

     # Axis limits for each fit parameter to plot
     ylims = [ (0, 1),
               (10**(-12), 10**(0)),
               (0.5, 0.9),
               (0, 200),
               (0, 200),
               (0, 8),
               (0, 0.4),
               (0, 0.2),
               (0, 0.6),
               (0, 0.25)]

     # Perform Mann-Whitney U test to compare R^2 distributions from the simulated
     # data sets
     print("p-values for Mann-Whitney U test between simulated R^2 distributions")
     ref_model_index = 0
     for m in range(len(colimitation_models.all_2D_trait_models)):
          if m == ref_model_index: continue
          Rsq1 = [fit_params_sims[s][ref_model_index][0] for s in range(num_sims) if not numpy.isnan(fit_params_sims[s][ref_model_index][0])]
          Rsq2 = [fit_params_sims[s][m][0] for s in range(num_sims) if not numpy.isnan(fit_params_sims[s][m][0])]
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
     
          # Plot fit parameters for simulated data
          params_sims = [[fit_params_sims[s][i][p] for s in range(num_sims) if not numpy.isnan(fit_params_sims[s][i][p])] for i in model_indices]
          boxes = axis.boxplot(params_sims, positions=model_indices, sym="", patch_artist=True, zorder=0, 
               medianprops={"color": "tab:purple", "linewidth": 0.5},
               boxprops={"facecolor": "tab:purple", "edgecolor": "tab:purple", "linewidth": 0.5},
               whiskerprops={"color": "tab:purple", "linewidth": 1.5},
               capprops={"color": "tab:purple", "linewidth": 1.5}, label="Data simulated from\n" + true_model_formatted + " model")

          # Mark true parameter values but skip first two plots (R^2 and Akaike weight)
          if p not in [0, 1]:
               axis.axhline(true_parameters[p - 2], linestyle="--", color="gray", label="True parameter value")
               axis.text(len(colimitation_models.all_2D_trait_models) - 1, 1.02*true_parameters[p - 2], r"\textbf{True value}", va="bottom", ha="right", fontsize=colimitation_plots.axis_label_size)

          # Show legend in first plot only
          if p == 0:
               axis.legend(loc=(0, 1.9), fontsize=colimitation_plots.axis_label_size)

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

     figure.savefig("Figure_S16.pdf", bbox_inches="tight")

################################################################################


if __name__ == '__main__':
     main()
