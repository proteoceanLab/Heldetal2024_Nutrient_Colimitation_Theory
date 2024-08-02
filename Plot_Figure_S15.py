import colimitation_data_analysis
import colimitation_plots
from matplotlib import pyplot
import pandas
import numpy
from scipy import stats


def main():

     # Number of replicate experiments
     num_replicates = 3

     # Number of R1 (glucose) and R2 (ammonium) values
     num_R1_values = 10
     num_R2_values = 7

     # Column labels for glucose and ammonium concentrations in data file
     R1_label = "Glucose (mM)"
     R2_label = "Ammonium (mM)"

     # Column labels for replicate rate measurements in data file
     rate_label_base = "Growth rate (per hour)"
     rate_labels = [f"{rate_label_base} rep {r + 1}" for r in range(num_replicates)]

     # Read R1, R2, and rate data into 2D meshes
     skiprows = 11
     sheet_name = "Table 1 Growth rates"
     R1_mesh, R2_mesh, rate_mesh = colimitation_data_analysis.Read2DScan("./Data/Dataset_S1.xlsx", R1_label, R2_label, rate_labels, num_R1_values, num_R2_values, sheet_name=sheet_name, skiprows=skiprows) 

     # Calculate linear dependence of standard deviation of rate with mean rate
     means = []     
     sds = []
     for i in range(len(rate_mesh)):
          for j in range(len(rate_mesh[i])):
               entries = [x for x in rate_mesh[i][j] if not isinstance(x, str)]
               if len(entries) > 0:
                    sds.append(numpy.std(entries))
                    means.append(numpy.mean(entries))
     result = stats.linregress(means, sds)

################################################################################

     # Set plot parameters
     pyplot.rcParams = colimitation_plots.SetrcParams(pyplot.rcParams)

     # Initialize figure
     figure = pyplot.figure(figsize=(4, 3))
     figure.subplots_adjust(hspace=0.05)
     
     # Plot means vs. standard deviations
     axis = figure.add_subplot(1, 1, 1)
     axis.set_xlabel("Mean growth rate across replicates (per hour)", fontsize=colimitation_plots.axis_label_size)
     axis.set_ylabel("Standard deviation of growth\nrate across replicates (per hour)", fontsize=colimitation_plots.axis_label_size)
     axis.set_xlim([0.1, 0.7])
     axis.set_ylim([0, 0.15])
     axis.scatter(means, sds, label="Data")
     axis.plot([min(means), max(means)], [result.intercept + result.slope*min(means), result.intercept + result.slope*max(means)], "--", color="red", 
          label="Linear regression" + "\n" + r"Slope $\approx$ " + str(round(result.slope, 3)) + "\n" + r"Intercept $\approx$ " + str(round(result.intercept, 3)) + "\n" + r"$p$-value $\approx$ " + str(round(result.pvalue, 3)))
     axis.legend(loc=1, fontsize=colimitation_plots.legend_label_size)

     figure.savefig("Figure_S15.pdf", bbox_inches="tight")

################################################################################


if __name__ == '__main__':
     main()
