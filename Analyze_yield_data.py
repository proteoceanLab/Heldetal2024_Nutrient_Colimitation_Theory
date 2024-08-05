import colimitation_data_analysis
import colimitation_models
import pandas


def main():
     """
     This script reads yield 2D scan data, fits it to all yield models, and 
     calculates virtual supplementation experiments.
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

     # Fit all replicates together and each replicate separately
     replicate_sets = [list(range(num_replicates))] + [[r] for r in range(num_replicates)]
     for reps in replicate_sets:
          fit_dataframe = colimitation_data_analysis.FitAllModels(R1_mesh, R2_mesh, yield_mesh, reps, colimitation_models.all_2D_trait_models)
          fit_dataframe.to_csv(f"./Data/Yield_data_fits_rep{''.join(str(r + 1) for r in reps)}.csv", index=False)

     # Iterate over replicates and calculate virtual supplementations
     for r in range(num_replicates):
          lim_coeffs_dataframe = colimitation_data_analysis.CalculateVirtualSupplementations(R1_mesh, R2_mesh, yield_mesh[:, :, r])
          lim_coeffs_dataframe.to_csv(f"./Data/Yield_data_lim_coeffs_rep{r + 1}.csv", index=False)

################################################################################

if __name__ == '__main__':
     main()