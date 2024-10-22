import colimitation_data_analysis
import colimitation_models
import pandas


def main():
     """
     This script generates bootstrapped yield data sets from the original 
     experimental data set and then fits them to all yield models.
     """

     # Number of replicate experiments
     num_replicates = 3

     # Number of R1 (glucose) and R2 (ammonium) values
     num_R1_values = 11
     num_R2_values = 7

     # Column labels for glucose and ammonium concentrations in data file
     R1_label = "Glucose (mM)"
     R2_label = "Ammonium (mM)"

     # Column labels for replicate yield measurements in data file
     yield_label_base = "Growth yield (OD 600 nm)"
     yield_labels = [f"{yield_label_base} rep {r + 1}" for r in range(num_replicates)]

     # Read R1, R2, and yield data into 2D meshes
     skiprows = 11
     sheet_name = "Table 1 Growth yields"
     R1_mesh, R2_mesh, yield_mesh = colimitation_data_analysis.Read2DScan("./Data/Dataset_S2.xlsx", R1_label, R2_label, yield_labels, num_R1_values, num_R2_values, sheet_name=sheet_name, skiprows=skiprows) 

     # Number of bootstrapped data sets
     num_bootstraps = 100

     # Seed for random number generation
     seed = 0

     # Generate bootstrapped data sets and write to files
     bootstraps_dataframe = colimitation_data_analysis.Generate2DScanBootstraps(R1_mesh, R2_mesh, yield_mesh, R1_label, R2_label, yield_label_base, num_bootstraps, seed=seed)
     bootstraps_dataframe.to_csv("./Data/Yield_bootstraps.csv", index=False)

     # Read back bootstrapped data from file
     yield_labels = [f"{yield_label_base} bootstrap {b + 1}" for b in range(num_bootstraps)]
     R1_mesh_bootstrap, R2_mesh_bootstrap, yield_mesh_bootstrap = colimitation_data_analysis.Read2DScan("./Data/Yield_bootstraps.csv", R1_label, R2_label, yield_labels, num_R1_values, num_R2_values) 

     # Set up combined dataframe for all bootstrap fits
     columns = ["Model name", "Num parameters", "Parameter names"] + [f"{c} bootstrap {b + 1}" for b in range(num_bootstraps) for c in ["R^2", "Akaike weight", "Parameter values"] ]
     fit_dataframe_combined = pandas.DataFrame([], columns=columns)

     # Iterate over each bootstrap
     for b in range(num_bootstraps):

          # Fit this bootstrap
          fit_dataframe = colimitation_data_analysis.FitAllModels(R1_mesh_bootstrap, R2_mesh_bootstrap, yield_mesh_bootstrap, [b], colimitation_models.all_2D_trait_models)

          # For the first bootstrap, add fixed columns to combined dataframe
          if b == 0:
               fit_dataframe_combined["Model name"] = fit_dataframe["Model name"]
               fit_dataframe_combined["Num parameters"] = fit_dataframe["Num parameters"]
               fit_dataframe_combined["Parameter names"] = fit_dataframe["Parameter names"]

          # Add columns specific to each bootstrap fit
          fit_dataframe_combined[f"R^2 bootstrap {b + 1}"] = fit_dataframe["R^2"]
          fit_dataframe_combined[f"Akaike weight bootstrap {b + 1}"] = fit_dataframe["Akaike weight"]
          fit_dataframe_combined[f"Parameter values bootstrap {b + 1}"] = fit_dataframe["Parameter values"]

     # Fit all bootstrap fit data to file
     fit_dataframe_combined.to_csv("./Data/Yield_bootstraps_fits.csv", index=False)

################################################################################

if __name__ == '__main__':
     main()