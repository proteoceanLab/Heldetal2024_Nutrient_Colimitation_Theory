import colimitation_data_analysis
import colimitation_models
import numpy
from scipy import stats
import pandas


def main():

     # Read experimental data

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
     R1_mesh, R2_mesh, yield_mesh = colimitation_data_analysis.Read2DScan("./Data/Dataset_S3_Rachana.xlsx", R1_label, R2_label, yield_labels, num_R1_values, num_R2_values, sheet_name=sheet_name, skiprows=skiprows) 

################################################################################

     # Determine linear regression of mean vs. standard deviation of yield
     # across replicates
     means = []     
     sds = []
     for i in range(len(yield_mesh)):
          for j in range(len(yield_mesh[i])):
               entries = [x for x in yield_mesh[i][j] if not isinstance(x, str)]
               if len(entries) > 0:
                    sds.append(numpy.std(entries))
                    means.append(numpy.mean(entries))
     result = stats.linregress(means, sds)

################################################################################

     # Simulate data sets from models with replicate noise based on actual data

     # Number of simulations per model
     num_sims = 10000

     # Seed for random number generation
     seed = 0

     # List of models to simulate
     models_to_simulate = ["liebig_blackman", "liebig_monod", "pat"]

     # Iterate over models
     for model_name in models_to_simulate:

          # Get parameters for that model from fits to the experimental data
          dataframe = pandas.read_csv("./Data/Yield_data_fits_rep123.csv", index_col=False)
          row = dataframe[dataframe["Model name"] == model_name]
          row_index = row.index[0]
          params = [float(x) for x in str(row["Parameter values"][row_index]).split(";")]

          # Simulate data for that model with noise based on experimental replicates
          model_function = lambda R1, R2: colimitation_models.CalcTraitForFit([R1, R2], params, model_name)
          error_function = lambda rate: rate*result.slope + result.intercept
          sims_dataframe = colimitation_data_analysis.Simulate2DScans(R1_mesh, R2_mesh, model_function, error_function, R1_label, R2_label, "Growth yield (OD 600 nm)", num_sims, seed)
          sims_dataframe.to_csv(f"./Data/Yield_sim_{model_name}.csv", index=False)
     
          # Read back simulated data from file
          yield_labels = [f"Growth yield (OD 600 nm) sim {s + 1}" for s in range(num_sims)]
          R1_mesh_sim, R2_mesh_sim, yield_mesh_sim  = colimitation_data_analysis.Read2DScan(f"./Data/Yield_sim_{model_name}.csv", R1_label, R2_label, yield_labels, num_R1_values, num_R2_values) 
     
          # Set up dataframe for fit parameters for all simulations
          columns = ["Model name", "Num parameters", "Parameter names"] + [f"{c} sim {s + 1}" for s in range(num_sims) for c in ["R^2", "Akaike weight", "Parameter values"] ]
          fit_dataframe_combined = pandas.DataFrame([], columns=columns)
     
          # Iterate over each simulation
          for s in range(num_sims):
     
               # Fit this simulation
               fit_dataframe = colimitation_data_analysis.FitAllModels(R1_mesh_sim, R2_mesh_sim, yield_mesh_sim, [s], colimitation_models.all_2D_trait_models)
     
               # For the first simulation, add fixed columns to combined dataframe
               if s == 0:
                    fit_dataframe_combined["Model name"] = fit_dataframe["Model name"]
                    fit_dataframe_combined["Num parameters"] = fit_dataframe["Num parameters"]
                    fit_dataframe_combined["Parameter names"] = fit_dataframe["Parameter names"]
     
               # Add columns specific to each simulation fit
               fit_dataframe_combined[f"R^2 sim {s + 1}"] = fit_dataframe["R^2"]
               fit_dataframe_combined[f"Akaike weight sim {s + 1}"] = fit_dataframe["Akaike weight"]
               fit_dataframe_combined[f"Parameter values sim {s + 1}"] = fit_dataframe["Parameter values"]
     
          # Write all simulation fit data to file
          fit_dataframe_combined.to_csv(f"./Data/Yield_sim_{model_name}_fits.csv", index=False)
     
          # Calculate virtual supplementations for noisy simulations
     
          # Range of thresholds of limitation coefficients for absolute colimitation test
          # List of lists of fractions of virtual supplementations with absolute colimitation above each threshold
          L_threshold_range = numpy.linspace(0, 0.2, 11)
          fracs_colim_sims = [[] for i in range(len(L_threshold_range))]
     
          # Iterate over simulations
          for s in range(num_sims):
     
               # Calculate limitation coefficients from virtual supplementations
               lim_coeffs_sim = colimitation_data_analysis.CalculateVirtualSupplementations(R1_mesh_sim, R2_mesh_sim, yield_mesh_sim[:, :, s])
     
               # Calculate fraction of virtual supplementations with absolute colimitation
               # above each threshold
               for i in range(len(L_threshold_range)):
                    frac = ((lim_coeffs_sim["L1"] > L_threshold_range[i]) & (lim_coeffs_sim["L2"] > L_threshold_range[i])).sum()/len(lim_coeffs_sim["Meff"])
                    fracs_colim_sims[i].append(frac)
     
          # Write distributions of fractions of absolute colimitation to file
          dataframe = pandas.DataFrame(numpy.array(fracs_colim_sims).T, columns=[f"Fraction of virtual supplementations with L1 and L2 > {t}" for t in L_threshold_range])
          dataframe.to_csv(f"./Data/Yield_sim_{model_name}_colim_fracs.csv", index=False)

################################################################################

     # Simulate data sets from models without noise --- only record limitation coefficients

     # Number of simulations per model
     num_sims = 1

     # Seed for random number generation
     seed = 0

     # List of models to simulate
     models_to_simulate = ["liebig_blackman", "liebig_monod", "pat"]

     # Iterate over models
     for model_name in models_to_simulate:

          # Get parameters for that model from fits to the experimental data
          dataframe = pandas.read_csv("./Data/Yield_data_fits_rep123.csv", index_col=False)
          row = dataframe[dataframe["Model name"] == model_name]
          row_index = row.index[0]
          params = [float(x) for x in str(row["Parameter values"][row_index]).split(";")]

          # Simulate data for that model with noise based on experimental replicates
          model_function = lambda R1, R2: colimitation_models.CalcTraitForFit([R1, R2], params, model_name)
          error_function = lambda rate: 0
          sims_dataframe = colimitation_data_analysis.Simulate2DScans(R1_mesh, R2_mesh, model_function, error_function, R1_label, R2_label, "Growth yield (OD 600 nm)", num_sims, seed)

          # Get meshes from simulated data
          R1_mesh_sim = numpy.array(sims_dataframe[R1_label]).reshape((num_R2_values, num_R1_values))
          R2_mesh_sim = numpy.array(sims_dataframe[R2_label]).reshape((num_R2_values, num_R1_values))
          yield_mesh_sim = numpy.array([numpy.array(sims_dataframe["Growth yield (OD 600 nm) sim 1"]).reshape((num_R2_values, num_R1_values))])
          yield_mesh_sim = numpy.transpose(yield_mesh_sim, (1, 2, 0))

          # Calculate limitation coefficients from virtual supplementations
          lim_coeffs_sim = colimitation_data_analysis.CalculateVirtualSupplementations(R1_mesh_sim, R2_mesh_sim, yield_mesh_sim[:, :, 0])
          lim_coeffs_sim.to_csv(f"./Data/Yield_sim_{model_name}_noiseless_lim_coeffs.csv", index=False)

################################################################################

if __name__ == '__main__':
     main()
