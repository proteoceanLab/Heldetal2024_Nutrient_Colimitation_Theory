import colimitation_data_analysis
import colimitation_models
import pandas


def main():
     """
     This script reads rate 2D scan data, fits it to all rate models, and 
     calculates virtual supplementation experiments.
     """

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

     # Fit all replicates together and each replicate separately
     replicate_sets = [list(range(num_replicates))] + [[r] for r in range(num_replicates)]
     for reps in replicate_sets:
          fit_dataframe = colimitation_data_analysis.FitAllModels(R1_mesh, R2_mesh, rate_mesh, reps, colimitation_models.all_2D_trait_models)
          fit_dataframe.to_csv(f"./Data/Rate_data_fits_rep{''.join(str(r + 1) for r in reps)}.csv", index=False)

     # Iterate over replicates and calculate virtual supplementations
     for r in range(num_replicates):

          # Get best-fit model and shift all R1 and R2 values by R1min and R2min for virtual supplementations
          dataframe = pandas.read_csv("./Data/Rate_data_fits_rep" + str(r + 1) + ".csv", index_col=False)
          row_index = dataframe["Akaike weight"].idxmax()
          best_model = dataframe.iloc[row_index]["Model name"]
          if "rmin" in best_model:
               params = [float(x) for x in str(dataframe.iloc[row_index]["Parameter values"]).split(";")]
               R1min, R2min = params[-2:]
          else:
               R1min, R2min = 0, 0
          R1_mesh_shifted = R1_mesh + R1min
          R2_mesh_shifted = R2_mesh + R2min

          # Calculate virtual supplementations
          lim_coeffs_dataframe = colimitation_data_analysis.CalculateVirtualSupplementations(R1_mesh_shifted, R2_mesh_shifted, rate_mesh[:, :, r])
          lim_coeffs_dataframe.to_csv(f"./Data/Rate_data_lim_coeffs_rep{r + 1}.csv", index=False)

################################################################################

if __name__ == '__main__':
     main()