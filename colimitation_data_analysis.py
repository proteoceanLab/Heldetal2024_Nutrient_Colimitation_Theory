"""
This module contains functions for analysis of growth rate and growth yield colimitation data (1D and 2D resource concentration scans).
"""


import numpy
import pandas
import random
from scipy import optimize
from matplotlib import pyplot
import colimitation_models


def Round(x, sigfigs):
     """
     x: Float
     sigfigs: Integer number of significant figures 

     return: String of the float x with the specified number of significant figures
     """

     return "%s" % float(("%." + str(sigfigs) + "g") % x)


def CalcRsq(x_values, data_values, function, params):
     """
     x_values: List of independent variable values
     data_values: List of dependent variable values
     function: A function that was fit to x_values and data_values (takes an x value and list of parameters as input)
     params: List of parameters for function

     return: R^2 of the fit
     """

     # Sum of squared residuals between data and the model
     sum_squared_residuals = sum((data_values - function(x_values, *params))**2)

     # Total sum of squares of between data and mean
     sum_squares = sum((data_values - numpy.mean(data_values))**2)

     return 1 - sum_squared_residuals/sum_squares


def CalcAICc(x_values, data_values, function, params):
     """
     x_values: List of independent variable values
     data_values: List of dependent variable values
     function: A function that was fit to x_values and data_values (takes an x value and list of parameters as input)
     params: List of parameters for function

     return: Corrected Akaike information criterion (AICc) of the fit

     This is an AIC calculation for a least-squares fit based on

     https://en.wikipedia.org/wiki/Akaike_information_criterion#Comparison_with_least_squares

     and Sec. 2.2 of Burnham and Anderson 2002.
     """

     # Sum of squared residuals between data and the model
     sum_squared_residuals = sum((data_values - function(x_values, *params))**2)

     # Initialize number of parameters and number of data points
     num_params = len(params)
     num_data_points = len(data_values)

     # Degrees of freedom
     dof = num_data_points - num_params

     # Calculate AIC
     aic = 2*num_params + num_data_points*numpy.log(sum_squared_residuals/(num_data_points - num_params))

     # This is a correction to AIC for AICc, accounting for small number of data points relative to number of parameters
     # (See Sec. 2.4 of Burnham and Anderson 2002)
     correction = 2*num_params*(num_params + 1)/(num_data_points - num_params - 1)
     aicc = aic - correction

     return aicc


def Read2DScan(file_name, R1_label, R2_label, trait_labels, num_R1_values, num_R2_values, sheet_name=0, skiprows=0):
     """
     file_name: Complete path to file name (must be Excel or CSV format)
     R1_label: Name of column for resource 1 concentrations
     R2_label: Name of column for resource 2 concentrations
     trait_labels: List of column names for traits to read (multiple replicates)
     sheet_name: Name of sheet to read for Excel files
     skip_rows: Number of rows to skip at the beginning

     return: 2D meshes of R1 values and R2 values (shapes are (num R2 values, num R1 values)) and 2D mesh of trait values (shape is (num R2 values, num R1 values, num replicates))
     """

     # Extract data from Excel or CSV file into a dataframe
     if ".xlsx" in file_name or ".xls" in file_name:
          data = pandas.read_excel(file_name, sheet_name=sheet_name, skiprows=skiprows, usecols=[R1_label, R2_label] + trait_labels)
     elif ".csv" in file_name:
          data = pandas.read_csv(file_name, names=[R1_label, R2_label] + trait_labels, header=0, skiprows=skiprows)
     else:
          raise ValueError("Unknown file format")

     # Set up meshes for R1 and R2 values; indices are R2 (row), R1 (column) to match wells on the plate
     R1_mesh = numpy.array(data[R1_label]).reshape((num_R2_values, num_R1_values))
     R2_mesh = numpy.array(data[R2_label]).reshape((num_R2_values, num_R1_values))

     # Set up mesh for trait values
     trait_mesh = numpy.array([numpy.array(data[trait_labels[i]]).reshape((num_R2_values, num_R1_values)) for i in range(len(trait_labels))])

     # Tranpose axes of mesh so that the indices are R2 (row), R1 (column), and replicate
     trait_mesh = numpy.transpose(trait_mesh, (1, 2, 0))

     return R1_mesh, R2_mesh, trait_mesh


def Convert2DScanDataForFit(R1_mesh, R2_mesh, trait_mesh, reps):
     """
     R1_mesh: 2D mesh of R1 values (shape is (num R2 values, num R1 values))
     R2_mesh: 2D mesh of R2 values (shape is (num R2 values, num R1 values))
     trait_mesh: 2D mesh of trait values (shape is (num R2 values, num R1 values, num replicates))
     reps: List of replicate indices to include in fits

     return: R1R2_values (list of R1 and R2 values with shape (2, num data points)) and trait_values (list of trait values with shape num data points)
     """

     # Initialize new flat arrays
     R1R2_values = []
     trait_values = []

     # Iterate over rows (R2 values)
     for i in range(len(trait_mesh)):

          # Iterate over columns (R1 values)
          for j in range(len(trait_mesh[i])):

               # Iterate over replicates to include
               for r in reps:

                    # Skip trait entries that aren't numbers (e.g., "not fit")
                    try:
                         trait_values.append(float(trait_mesh[i][j][r]))
                         R1R2_values.append([R1_mesh[i][j], R2_mesh[i][j]])
                    except ValueError:
                         pass

     # Convert to Numpy arrays and transpose for fitting
     R1R2_values = numpy.array(R1R2_values).T
     trait_values = numpy.array(trait_values)

     return R1R2_values, trait_values


def FitData(R_values, trait_values, model):
     """
     R_values: List of resource concentration values with dimensions (num_resources, num_data_points)
     trait_values: List of trait (either growth rate or biomass yield) values
     model: name of model to fit

     return: fitted parameter values, standard deviation of those fitted parameter values, R^2, and AICc
     """

     # Initialize guesses for parameter values

     # The first parameter in all models is the max trait so guess the maximum measured value
     param_guesses = [max(trait_values)]

     # Set guesses for other parameters depending on the model
     if len(R_values) == 1:

          # Half the resource concentration range: use this to set initial guess
          # for scale of resource concentration
          half_R1_range = 0.5*(max(R_values[0]) + min(R_values[0]))

          # The next parameters are always the trait to resource coefficients,
          # so guess the ratios of the max trait to half resource concentration range
          param_guesses += [max(trait_values)/half_R1_range]
          # For the Hill model, guess a Hill coefficient of 1 (equivalent to Monod)
          if "hill" in model:
               param_guesses += [1]
          # For rmin models, guess Rmin is zero
          if "rmin" in model:
               param_guesses += [0]

     elif len(R_values) == 2:

          # Half the resource concentration range: use this to set initial guess
          # for scale of resource concentration
          half_R1_range = 0.5*(max(R_values[0]) + min(R_values[0]))
          half_R2_range = 0.5*(max(R_values[1]) + min(R_values[1]))

          # The next parameters are always the trait to resource coefficients,
          # so guess the ratios of the max trait to half resource concentration range
          param_guesses += [max(trait_values)/half_R1_range, max(trait_values)/half_R2_range]
          # For the Hill model, guess a Hill coefficient of 1 (equivalent to Monod)
          if "hill" in model:
               param_guesses += [1, 1]
          # For the generalized-additive model, guess a large q value that makes it Liebig-like
          if "generalized_additive" in model:
               param_guesses += [colimitation_models.LIEBIG_LIKE_EXPONENT]
          # For rmin models, guess Rmin is zero
          if "rmin" in model:
               param_guesses += [0, 0]

     else:
          raise ValueError("Too many resources for fitting")

     # Set all lower bounds to zero, all upper bounds to infinity
     # Note that for the generalized-additive model of yield, this requires it to be
     # parameterized such that q > 0 means essential resources
     num_params = len(param_guesses)
     bounds = (num_params*[0], num_params*[numpy.inf])

     # Try the fit
     try:
          # Define wrapper function
          function = lambda Rs, *params: colimitation_models.CalcTraitForFit(Rs, params, model)

          # Run fit optimization
          params, covs = optimize.curve_fit(function, R_values, trait_values, p0=param_guesses, bounds=bounds)

          # Calculate standard deviations for fitted parameters
          sds = numpy.sqrt(numpy.diag(covs))

          # Calculate R^2 values for fits
          Rsq = CalcRsq(R_values, trait_values, function, params)

          # Calculate AIC for fit
          aicc = CalcAICc(R_values, trait_values, function, params)
     except RuntimeError:
          params = num_params*["NA"]
          sds = num_params*["NA"]
          Rsq = "NA"
          aicc = "NA"

     return params, sds, Rsq, aicc


def PlotFit(file_name, R_values, trait_values, R_labels, trait_label, model_function, params, data_only=False, interactive=False):
     """
     file_name: Root of file name to save plot PDF
     R_values: dimensions (number of resources, number of data points)
     trait_values: dimensions (number of data points)
     Rs_labels: List of axis labels for resources
     trait_label: Axis label for trait 
     model_function: Fitted function to plot against data (call signature must be model_function(Rs, params))
     params: Parameters to pass to model_function
     data_only: If true, plot only the data and not the model
     interactive: If true, display plot in interactive mode
     """

     figure = pyplot.figure(figsize=(4, 3))

     # Plot 1D for a single resource
     if len(R_values) == 1:
          axis = figure.add_subplot(1, 1, 1)
          axis.set_xlabel(R_labels[0])
          axis.set_ylabel(trait_label)

          # Plot data
          axis.scatter(R_values[0], trait_values, color="tab:red")

          # Plot curve of fitted model
          if not data_only:
               R1_range = numpy.linspace(0, max(R_values[0]), 50)
               trait_range = model_function((R1_range), params)
               axis.plot(R1_range, trait_range, "-", color="tab:blue")

     # Plot 2D for two resources
     elif len(R_values) == 2:

          # Plot dummy axis first to make sure 3D plot isn't cut off
          axis = figure.add_subplot(1, 1, 1)
          axis.axis("off")
          axis.text(-0.1, 0.5, ".", alpha=0)
          axis.text(0.5, -0.1, ".", alpha=0)

          # Plot actual 3D axis
          axis = figure.add_subplot(1, 1, 1, projection="3d")
          axis.set_xlabel(R_labels[0])
          axis.set_ylabel(R_labels[1])
          axis.set_zlabel(trait_label)
          axis.view_init(30, -135)
     
          # Plot data
          axis.scatter(R_values[0], R_values[1], trait_values, color="tab:red")
     
          # Plot surface of fitted model
          if not data_only:
               R1_range = numpy.linspace(0, max(R_values[0]), 50)
               R2_range = numpy.linspace(0, max(R_values[1]), 50)
               R1_mesh, R2_mesh = numpy.meshgrid(R1_range, R2_range)
               trait_mesh = model_function((R1_mesh, R2_mesh), params)
               axis.plot_surface(R1_mesh, R2_mesh, trait_mesh)
     else:
          raise ValueError("Too many resources for plotting")

     figure.savefig(file_name + ".pdf", bbox_inches="tight")
     if interactive:
          pyplot.show()
     pyplot.close()


def FitAllModels(R1_mesh, R2_mesh, trait_mesh, reps, list_of_models):
     """
     R1_mesh: 2D mesh of R1 values (shape is (num R2 values, num R1 values))
     R2_mesh: 2D mesh of R2 values (shape is (num R2 values, num R1 values))
     trait_mesh: 2D mesh of trait values (shape is (num R2 values, num R1 values, num replicates))
     reps: List of replicate indices to include in fits
     list_of_models: List of model names to fit

     return: dataframe with all fit data
     """

     # Convert 2D scan meshes to lists for fit
     R1R2_values, trait_values = Convert2DScanDataForFit(R1_mesh, R2_mesh, trait_mesh, reps=reps)

     # Table of fitting data across all models
     fit_data = []

     # Iterate over models
     for model in list_of_models:

          # Fit the model
          params, sds, Rsq, aicc = FitData(R1R2_values, trait_values, model)

          # Set string of parameter names
          param_names = "zmax;c1;c2"
          if "hill" in model:
               param_names += ";n1;n2"
          if "generalized_additive" in model:
               param_names += ";q"
          if "rmin" in model:
               param_names += ";R1min;R2min"

          # Set string of parameter values
          param_string = ";".join(str(x) for x in params)

          fit_data.append([model, len(params), Rsq, aicc, param_names, param_string])

     # Convert AICc values to normalized Akaike weights
     aiccs = [row[3] for row in fit_data]
     min_aicc = min([a for a in aiccs if not isinstance(a, str)])
     akaike_weights = [numpy.exp(-0.5*(a - min_aicc)) if not isinstance(a, str) else 0 for a in aiccs]
     akaike_weights = akaike_weights/sum(akaike_weights)
     fit_data = [fit_data[i][:3] + [akaike_weights[i]] + fit_data[i][4:] for i in range(len(fit_data))]

     # Create dataframe with all the fit data
     fit_dataframe = pandas.DataFrame(fit_data, columns=["Model name", "Num parameters", "R^2", "Akaike weight", "Parameter names", "Parameter values"])

     return fit_dataframe


def CalculateVirtualSupplementations(R1_mesh, R2_mesh, trait_mesh):
     """
     R1_mesh: 2D mesh of R1 values with shape (num R2 values, num R1 values)
     R2_mesh: 2D mesh of R2 values with shape (num R2 values, num R1 values)
     trait_mesh: 2D mesh of trait values with shape (num R2 values, num R1 values)

     return: List of (L1, L2, Meff) combinations for all virtual supplementations
     """

     # List of matched L1,L2,Meff values over all virtual supplementation experiments
     lim_coeff_list = []

     # Iterate over background concentrations
     for i_back in range(len(trait_mesh)):
          for j_back in range(len(trait_mesh[i_back])):

               # Background concentrations and trait on which to perform supplementation
               R1_back = R1_mesh[i_back][j_back]
               R2_back = R2_mesh[i_back][j_back]
               trait_back = trait_mesh[i_back][j_back]

               # If background trait is close to zero or NA, skip because limitation coefficient is undefined
               if isinstance(trait_back, str): continue
               if trait_back < 1e-8: continue
               if R1_back < 1e-8 and R2_back < 1e-8: continue

               # Supplement this background concentration to all higher R2 values
               for i_supp in range(i_back + 1, len(trait_mesh)):

                    # Supplemented R2 concentration and trait 
                    R2_supp = R2_mesh[i_supp][j_back]
                    trait_supp = trait_mesh[i_supp][j_back]

                    if isinstance(trait_supp, str): continue

                    # Calculate limitation coefficient from R2 supplement
                    L2 = ( (trait_supp - trait_back)/(R2_supp - R2_back) )*( R2_back/trait_back )

                    # Supplement this background concentration to all higher R1 values and compare to R2 supplement
                    for j_supp in range(j_back + 1, len(trait_mesh[i_back])):

                         # Supplemented R1 concentration and trait
                         R1_supp = R1_mesh[i_back][j_supp]
                         trait_supp = trait_mesh[i_back][j_supp]

                         if isinstance(trait_supp, str): continue

                         # Calculate limitation coefficient from R1 supplement
                         L1 = ( (trait_supp - trait_back)/(R1_supp - R1_back) )*( R1_back/trait_back )

                         # Calculate Meff for this pair of supplements
                         if abs(L1) > 1e-8 or abs(L2) > 1e-8:
                              Meff = (L1 + L2)/max([abs(L1), abs(L2)])
                         else:
                              Meff = "NA"

                         # Add L1, L2, Meff for this supplement combination to list
                         lim_coeff_list.append((R1_back, R2_back, trait_back, R1_supp, R2_supp, L1, L2, Meff))

     # Create dataframe with limitation coefficient data
     dataframe = pandas.DataFrame(lim_coeff_list, columns=["R1_back", "R2_back", "trait_back", "R1_supp", "R2_supp", "L1", "L2", "Meff"])

     return dataframe


def Generate2DScanBootstraps(R1_mesh, R2_mesh, trait_mesh, R1_label, R2_label, trait_label, num_bootstraps, seed=None):
     """
     R1_mesh: shape (num R2 values, num R1 values)
     R2_mesh: shape (num R2 values, num R1 values)
     trait_mesh: shape (num R2 values, num R1 values, num replicates)
     R1_label: Name of column for resource 1 concentrations
     R2_label: Name of column for resource 2 concentrations
     trait_label: Name of column for trait
     num_bootstraps: Number of bootstrapped data sets to generate
     seed: Seed for random number generator
     """

     # Seed random number generator
     random.seed(seed)

     # List of bootstrapped data
     data = []

     # Iterate over concentration points and resample replicate traits
     for i in range(len(trait_mesh)):
          for j in range(len(trait_mesh[i])):
               resampled_data = random.choices(trait_mesh[i][j], k=num_bootstraps)
               data.append([R1_mesh[i][j], R2_mesh[i][j], *resampled_data])

     # Create dataframe with all the bootstrapped data
     dataframe = pandas.DataFrame(data, columns=[R1_label, R2_label] + [f"{trait_label} bootstrap {b + 1}" for b in range(num_bootstraps)])

     return dataframe


def Simulate2DScans(R1_mesh, R2_mesh, model_function, error_function, R1_label, R2_label, trait_label, num_sims, seed=None):
     """
     R1_mesh: shape (num R2 values, num R1 values)
     R2_mesh: shape (num R2 values, num R1 values)
     model_function: function that converts R1, R2 into a trait (rate or yield)
     error_function: function that takes the mean trait value and determines the
          corresponding standard deviation for a Gaussian model
     R1_label: Name of column for resource 1 concentrations
     R2_label: Name of column for resource 2 concentrations
     trait_label: Name of column for trait
     num_sims: Number of simulations
     seed: Seed for random number generator
     """

     # Seed random number generator
     random.seed(seed)

     # List of simulated data
     data = []

     # Iterate over concentration points and simulate trait values
     for i in range(len(R1_mesh)):
          for j in range(len(R1_mesh[i])):

               # Calculate mean and standard deviation of trait at this point
               trait_mean = model_function(R1_mesh[i][j], R2_mesh[i][j])
               trait_sd = error_function(trait_mean)

               # Generate simulated trait values
               simulated_data = [random.gauss(trait_mean, trait_sd) for s in range(num_sims)]
               data.append([R1_mesh[i][j], R2_mesh[i][j], *simulated_data])

     # Create dataframe with all the bootstrapped data
     dataframe = pandas.DataFrame(data, columns=[R1_label, R2_label] + [f"{trait_label} sim {s + 1}" for s in range(num_sims)])

     return dataframe