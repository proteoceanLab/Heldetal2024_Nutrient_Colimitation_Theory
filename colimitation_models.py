"""
This module contains functions for models related to colimitation of multiple resources.
"""

import numpy
from scipy import integrate
from scipy import interpolate
from fast_helper_functions import CalcPATSum
from fast_helper_functions import CalcRateLimitationPAT
import warnings


# Exponent to use in generalized-mean model that approximates Liebig regime
# Note that this assumes the generalized-mean model is parameterized such that 
# q > 0 means essential resources
LIEBIG_LIKE_EXPONENT = 20


all_2D_trait_models = [  "liebig_monod", 
                         "liebig_monod_rmin",                   
                         "liebig_hill",                         
                         "liebig_hill_rmin",                    
                         "liebig_blackman",                     
                         "liebig_blackman_rmin",                
                         "liebig_bertalanffy",                  
                         "liebig_bertalanffy_rmin",             
                         "multiplicative_monod",                 
                         "multiplicative_monod_rmin",           
                         "multiplicative_hill",                  
                         "multiplicative_hill_rmin",             
                         "multiplicative_blackman",              
                         "multiplicative_blackman_rmin",         
                         "multiplicative_bertalanffy",           
                         "multiplicative_bertalanffy_rmin",      
                         "pat",      
                         "pat_rmin",
                         "additive", 
                         "additive_rmin",
                         "mankad_bungay",                        
                         "mankad_bungay_rmin",
                         "generalized_mean",
                         "generalized_mean_rmin",
                         "saito_sub",
                         "saito_sub_rmin",                       
                         "mean_monod_sub",                      
                         "mean_monod_sub_rmin",                 
                         "chem_depend",  
                         "chem_depend_rmin" 
]

all_2D_trait_models_formatted = [" ".join(word.capitalize() for word in model.split("_")).replace("Pat", "PAT") for model in all_2D_trait_models]


def CalcTraitODE(Rs, ss, zmax, model):
     """
     Rs: List of resource concentrations
     ss: List of stoichiometry coefficients
     zmax: Maximum trait value
     model: Name of trait model to use

     return: Trait value with the given parameters and model

     This version of the trait model calculation is optimized for use in ODEs
     with an arbitrary number of resources.
     """

     # Standardized parameters
     Rs = numpy.array(Rs)
     ss = numpy.array(ss)
     model = model.lower()

     # Set dimensionless resource concentrations
     rs = ss*Rs/zmax

     # Pre-factor to ensure traits are always zero (for essential independent 
     # resource models) if any resource concentration drops below zero
     no_zeros = int(numpy.all(rs > 0))

     # Calculate Liebig growth rate
     if model == "liebig_monod":
          return no_zeros*zmax*numpy.min(rs/(rs + 1))

     # Calculate smooth (approximate) Liebig growth rate
     elif model == "liebig_smooth":
          fluxes = zmax*rs/(rs + 1)
          return no_zeros*sum(fluxes**(-LIEBIG_LIKE_EXPONENT))**(1/(-LIEBIG_LIKE_EXPONENT))

     # Calculate the Poisson arrival time (PAT) growth rate
     elif model == "pat":

          # Calculate the complicated sum in the denominator using a separate, compiled Cython function
          total = CalcPATSum(rs)

          return no_zeros*zmax/(1 + total)

     # Calculate the additive growth rate (formerly known as the sequential-uptake model)
     elif model == "additive":
          return no_zeros*zmax/(1 + sum(1/rs))

     # Calculate the multiplicative growth rate
     elif model == "multiplicative_monod":
          return no_zeros*zmax*numpy.prod(rs/(rs + 1))

     # Calculate the Saito et al. model of substitutable resources
     elif model == "saito_sub":
          return zmax*sum(rs)/(1 + sum(rs))

     # Calculate the mean Monod model of substitutable resources
     elif model == "mean_monod_sub":
          return zmax*numpy.mean(rs/(rs + 1))

     # Calculate the Saito et al. model of chemically-dependent resources
     elif model == "chem_depend":
          if len(Rs) != 2:
               raise ValueError("Invalid number of resources for the chemically-dependent model")
          return zmax*Rs[0]/(Rs[0] + Ks[0]*(Rs[1] + Ks[1])/Rs[1])

     # If model is unknown, raise an error
     else:
          raise ValueError("Unknown growth rate model")


def CalcTraitForFit(Rs, params, model):
     """
     Rs: List of resource concentrations
     params: List of parameter values (specific to the model)
     model: Name of trait model to use

     return: Trait value under the given model and parameters

     This version of the trait model calculation is specialized to one or two 
     resources only and designed to be used with a fitting algorithm.
     """

     # Standardize model name to lowercase
     model = model.lower()

     # Small number to test for avoiding divisions by zero
     epsilon = 1e-10

     # Ignore runtime warnings for this function since they often result from
     # divide by zeros
     with warnings.catch_warnings():
          warnings.filterwarnings("ignore", message="divide by zero encountered in reciprocal", category=RuntimeWarning)
          warnings.filterwarnings("ignore", message="divide by zero encountered in power", category=RuntimeWarning)
          warnings.filterwarnings("ignore", message="divide by zero encountered in scalar power", category=RuntimeWarning)
          warnings.filterwarnings("ignore", message="overflow encountered in power", category=RuntimeWarning)
          warnings.filterwarnings("ignore", message="overflow encountered in add", category=RuntimeWarning)
     
          # If fitting 1D data:
          if len(Rs) == 1:
               R, = Rs
     
               if model == "monod":
                    zmax, s = params
                    return zmax*s*R/(s*R + zmax)
     
               elif model == "monod_rmin":
                    zmax, s, Rmin = params
                    return CalcTraitForFit([R + Rmin], [zmax, s], "monod")
     
               elif model == "hill":
                    zmax, s, n = params
                    return zmax*(s*R)**n/((s*R)**n + zmax**n)
     
               elif model == "hill_rmin":
                    zmax, s, n, Rmin = params
                    return CalcTraitForFit([R + Rmin], [zmax, s, n], "hill")
     
               elif model == "blackman":
                    zmax, s = params
                    return numpy.minimum(s*R, zmax)
     
               elif model == "blackman_rmin":
                    zmax, s, Rmin = params
                    return CalcTraitForFit([R + Rmin], [zmax, s], "blackman")
     
               elif model == "bertalanffy":
                    zmax, s = params
                    return zmax*(1 - 2**(-s*R/zmax))
     
               elif model == "bertalanffy_rmin":
                    zmax, s, Rmin = params
                    return CalcTraitForFit([R + Rmin], [zmax, s], "bertalanffy")
     
               else:
                    raise ValueError("Unknown 1D trait model")
     
          # If fitting 2D data:
          elif len(Rs) == 2:
               R1, R2 = Rs
     
               # Liebig models
               if model == "liebig_monod":
                    zmax, s1, s2 = params
                    z1 = CalcTraitForFit([R1], [zmax, s1], "monod")
                    z2 = CalcTraitForFit([R2], [zmax, s2], "monod")
                    return numpy.minimum(z1, z2)
     
               elif model == "liebig_monod_rmin":
                    zmax, s1, s2, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2], "liebig_monod")
     
               elif model == "liebig_hill":
                    zmax, s1, s2, n1, n2 = params
                    z1 = CalcTraitForFit([R1], [zmax, s1, n1], "hill")
                    z2 = CalcTraitForFit([R2], [zmax, s2, n2], "hill")
                    return numpy.minimum(z1, z2)
     
               elif model == "liebig_hill_rmin":
                    zmax, s1, s2, n1, n2, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2, n1, n2], "liebig_hill")
     
               elif model == "liebig_blackman":
                    zmax, s1, s2 = params
                    z1 = CalcTraitForFit([R1], [zmax, s1], "blackman")
                    z2 = CalcTraitForFit([R2], [zmax, s2], "blackman")
                    return numpy.minimum(z1, z2)
     
               elif model == "liebig_blackman_rmin":
                    zmax, s1, s2, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2], "liebig_blackman")
     
               elif model == "liebig_bertalanffy":
                    zmax, s1, s2 = params
                    z1 = CalcTraitForFit([R1], [zmax, s1], "bertalanffy")
                    z2 = CalcTraitForFit([R2], [zmax, s2], "bertalanffy")
                    return numpy.minimum(z1, z2)
     
               elif model == "liebig_bertalanffy_rmin":
                    zmax, s1, s2, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2], "liebig_bertalanffy")
     
               # Multiplicative models
               elif model == "multiplicative_monod":
                    zmax, s1, s2 = params
                    z1 = CalcTraitForFit([R1], [zmax, s1], "monod")
                    z2 = CalcTraitForFit([R2], [zmax, s2], "monod")
                    return z1*z2/zmax
     
               elif model == "multiplicative_monod_rmin":
                    zmax, s1, s2, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2], "multiplicative_monod")
     
               elif model == "multiplicative_hill":
                    zmax, s1, s2, n1, n2 = params
                    z1 = CalcTraitForFit([R1], [zmax, s1, n1], "hill")
                    z2 = CalcTraitForFit([R2], [zmax, s2, n2], "hill")
                    return z1*z2/zmax
     
               elif model == "multiplicative_hill_rmin":
                    zmax, s1, s2, n1, n2, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2, n1, n2], "multiplicative_hill")
     
               elif model == "multiplicative_blackman":
                    zmax, s1, s2 = params
                    z1 = CalcTraitForFit([R1], [zmax, s1], "blackman")
                    z2 = CalcTraitForFit([R2], [zmax, s2], "blackman")
                    return z1*z2/zmax
     
               elif model == "multiplicative_blackman_rmin":
                    zmax, s1, s2, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2], "multiplicative_blackman")
     
               elif model == "multiplicative_bertalanffy":
                    zmax, s1, s2 = params
                    z1 = CalcTraitForFit([R1], [zmax, s1], "bertalanffy")
                    z2 = CalcTraitForFit([R2], [zmax, s2], "bertalanffy")
                    return z1*z2/zmax
     
               elif model == "multiplicative_bertalanffy_rmin":
                    zmax, s1, s2, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2], "multiplicative_bertalanffy")
     
               # PAT, additive, and Mankad-Bungay models
               elif model == "pat":
                    zmax, s1, s2 = params
                    denominator = numpy.array((s1*R1)*(s2*R2)*(s1*R1 + s2*R2) + zmax*((s1*R1)**2 + (s1*R1)*(s2*R2) + (s2*R2)**2))
                    # Mask denominator values below a small number to avoid division by zero; this doesn't change function values because these points must be zero anyway
                    denominator[denominator < epsilon] = epsilon
                    return zmax*(s1*R1)*(s2*R2)*(s1*R1 + s2*R2)/denominator
     
               elif model == "pat_rmin":
                    zmax, s1, s2, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2], "pat")
     
               elif model == "additive":
                    zmax, s1, s2 = params
                    denominator = numpy.array((s1*R1)*(s2*R2) + zmax*(s1*R1 + s2*R2))
                    # Mask denominator values below a small number to avoid division by zero; this doesn't change function values because these points must be zero anyway
                    denominator[denominator < epsilon] = epsilon
                    return zmax*(s1*R1)*(s2*R2)/denominator
     
               elif model == "additive_rmin":
                    zmax, s1, s2, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2], "additive")
     
               elif model == "mankad_bungay":
                    zmax, s1, s2 = params
                    denominator = numpy.array(s1*R1 + s2*R2)
                    # Mask denominator values below a small number to avoid division by zero; this doesn't change function values because these points must be zero anyway
                    denominator[denominator < epsilon] = epsilon
                    return zmax*((s1*R1)*(s2*R2)/denominator)*(1/(s1*R1 + zmax) + 1/(s2*R2 + zmax))
     
               elif model == "mankad_bungay_rmin":
                    zmax, s1, s2, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2], "mankad_bungay")
     
               elif model == "generalized_mean":
                    zmax, s1, s2, q = params
                    return ((s1*R1)**(-q) + (s2*R2)**(-q) + zmax**(-q))**(1/(-q))
     
               elif model == "generalized_mean_rmin":
                    zmax, s1, s2, q, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2, q], "generalized_mean")
     
               # Models for substitutable and chemically-dependent resources
               elif model == "saito_sub":
                    zmax, s1, s2 = params
                    return zmax*(s1*R1 + s2*R2)/(s1*R1 + s2*R2 + zmax)
     
               elif model == "saito_sub_rmin":
                    zmax, s1, s2, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2], "saito_sub")
     
               elif model == "mean_monod_sub":
                    zmax, s1, s2 = params
                    return zmax*0.5*(s1*R1/(s1*R1 + zmax) + s2*R2/(s2*R2 + zmax))
     
               elif model == "mean_monod_sub_rmin":
                    zmax, s1, s2, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2], "mean_monod_sub")
     
               elif model == "chem_depend":
                    zmax, s1, s2 = params
                    return zmax*(s1*R1)*(s2*R2)/((s1*R1)*(s2*R2) + zmax*(s2*R2 + zmax))
     
               elif model == "chem_depend_rmin":
                    zmax, s1, s2, R1min, R2min = params
                    return CalcTraitForFit([R1 + R1min, R2 + R2min], [zmax, s1, s2], "chem_depend")
     
               else:
                    raise ValueError("Unknown 2D trait model: " + model)
     
          # If fitting more than 2 resource, raise an error for now
          else:
               raise ValueError("Too many resources for trait model")


def CalcLimCoeff(resource_index, Rs, ss, zmax, model):
     """
     resource_index: Index of resource for which to calculate limitation coefficient
     Rs: List of resource concentrations
     ss: List of stoichiometry coefficients
     zmax: Maximum trait value
     model: Name of trait model to use

     return: Limitation coefficient with the given parameters and model
     """

     # Standardized parameters
     Rs = numpy.array(Rs)
     ss = numpy.array(ss)
     model = model.lower()

     # Set dimensionless resource concentrations
     rs = ss*Rs/zmax

     # If any concentrations are negative (from integrating ODEs), then assume
     # the real concentrations are zero.  For all models except multiplicative,
     # this means the limitation coefficients for those resouces are 1, and 0 
     # for all nonzero resources.  For the multiplicative model, the limitation 
     # coefficients for zero resources are 1 but the limitation coefficients for
     # nonzero resources are as usual.

     # If this resource is zero, then its limitation is 1 for any model
     if rs[resource_index] < 0:
          return 1
     # If this resource is nonzero, but if another resource is zero, then its
     # limitation depends on the model
     elif any(rs < 0):
          if model in ["liebig_monod", "liebig_smooth", "pat", "additive"]:
               return 0
     
     # Calculate Liebig limitation coefficient
     if model == "liebig_monod":
          return int(resource_index == numpy.argmin(rs/(rs + 1)))/(rs[resource_index] + 1)

     # Calculate smooth (approximate) Liebig limitation coefficient
     elif model == "liebig_smooth":
          fluxes = rs/(rs + 1)
          return (1/(rs[resource_index] + 1))*fluxes[resource_index]**(-LIEBIG_LIKE_EXPONENT)/sum(fluxes**(-LIEBIG_LIKE_EXPONENT))

     # Calculate the Poisson arrival time (PAT) limitation coefficient
     elif model == "pat":
          # Calculate the complicated things using a separate, compiled Cython function
          return CalcRateLimitationPAT(resource_index, rs)

     # Calculate the additive limitation coefficient
     elif model == "additive":
          return (1/rs[resource_index])/(1 + sum(1/rs))

     # Calculate the multiplicative limitation coefficient
     elif model == "multiplicative_monod":
          return 1/(rs[resource_index] + 1)

     # Calculate the Saito et al. model of substitutable resources
     elif model == "saito_sub":
          return rs[resource_index]/sum(rs)/(1 + sum(rs))

     # Calculate the mean Monod model of substitutable resources
     elif model == "mean_monod_sub":
          return (rs[resource_index]/(rs[resource_index] + 1)**2/len(rs))/numpy.mean(rs/(rs + 1))

     # Calculate the Saito et al. model of chemically-dependent resources
     elif model == "chem_depend":
          if len(rs) != 2:
               raise ValueError("Invalid number of resources for the chemically-dependent model")
          r1, r2 = rs
          if resource_index == 0:
               return (r2 + 1)/(r1*r2 + r2 + 1)
          else:
               return 1/(r1*r2 + r2 + 1)

     # If model is unknown, raise an error
     else:
          raise ValueError("Unknown growth rate model")


def CalcLimCoeffGM(resource_index, Rs, ss, q):
     """
     resource_index: Index of resource for which to calculate limitation coefficient
     Rs: List of resource concentrations
     ss: List of stoichiometry coefficients
     q: Interaction parameter

     return: Limitation coefficient specifically for generalized mean model 
          (no implicit factors)
     """

     Rs = numpy.array(Rs)
     ss = numpy.array(ss)

     rs = ss*Rs

     return (rs[resource_index])**(-q)/sum(rs**(-q))


def CalcMeff(Rs, ss, zmax, model):
     """
     Rs: List of resource concentrations
     ss: List of stoichiometry coefficients
     zmax: Maximum trait value
     model: Name of trait model to use

     return: Number of effectively limiting resources
     """

     # Standardized parameters
     model = model.lower()

     # Calculate list of all limitation coefficients for explicit resources
     Lis = [CalcLimCoeff(i, Rs, ss, zmax, model) for i in range(len(Rs))]

     # Add limitation coefficient for implicit resource
     Lis.append(1 - sum(Lis))

     return 1/numpy.max(numpy.array(Lis))


def CalcMeffGM(Rs, ss, q):
     """
     Rs: List of resource concentrations
     ss: List of stoichiometry coefficients
     q: Interaction parameter

     return: Number of effectively limiting resources specifically for generalized mean model
          (no implicit factors)
     """

     # Calculate list of all limitation coefficients for explicit resources
     Lis = [CalcLimCoeffGM(i, Rs, ss, q) for i in range(len(Rs))]

     return 1/numpy.max(numpy.array(Lis))


def CalcSerialTransferSteadyState(Rsources, Ys, D):
     """
     Rsources: List of source concentrations for all resources
     Ys: List of intrinsic yields (same length as Rsources)
     D: Dilution factor (> 1)

     return: Initial concentrations of biomass and all resources in steady-state
          of serial transfers
     """

     # Standardized parameters
     Rsources = numpy.array(Rsources)
     Ys = numpy.array(Ys)

     # Identify the yield-limiting resource
     yield_limiting_resource = numpy.argmin(Rsources*Ys)

     # Calculate the biomass at the beginning of each cycle
     N0 = Rsources[yield_limiting_resource]*Ys[yield_limiting_resource]/(D - 1)
     N0s = numpy.array([N0])

     # Calculate the resource concentrations at the beginning of each cycle
     R0s = (D - Rsources[yield_limiting_resource]*Ys[yield_limiting_resource]/(Rsources*Ys))*Rsources/(D - 1)

     return N0s, R0s


def CalcRStarAdditiveModel(gmax, K1, K2, d, R1initial, R2initial, Y1, Y2):
     """
     gmax: Maximum growth rate in additive model
     K1, K2: Half-saturation concentrations in additive model
     d: Death rate or minimum growth rate
     R1initial, R2initial: Initial concentrations of resources
     Y1, Y2: Yields of resources

     return: Rstar values for resource 1 and resource 2
     """

     # Define dimensionless quantities
     gamma = gmax/d
     r1initial = R1initial/K1
     r2initial = R2initial/K2
     y1 = K1*Y1
     y2 = K2*Y2
     D = y1**2*(1 + r1initial*(r1initial*(gamma - 1) - 2)*(gamma - 1)) + y2**2*(1 + r2initial*(r2initial*(gamma - 1) - 2)*(gamma - 1)) + 2*y1*y2*(1 - r1initial - r2initial - r1initial*r2initial + (r1initial + r2initial + 2*r1initial*r2initial)*gamma - r1initial*r2initial*gamma**2)

     # Calculate Rstar values
     R1star = (y1 + y2 + (gamma - 1)*(y1*r1initial - y2*r2initial) + numpy.sqrt(D))*K1/(2*y1*(gamma - 1))
     R2star = (y1 + y2 + (gamma - 1)*(y2*r2initial - y1*r1initial) + numpy.sqrt(D))*K2/(2*y2*(gamma - 1))

     return R1star, R2star


def CalcYieldChemostatAdditiveModel(gmax, K1, K2, d, R1initial, R2initial, Y1, Y2):
     """
     gmax: Maximum growth rate in additive model
     K1, K2: Half-saturation concentrations in additive model
     d: Death rate or minimum growth rate
     R1initial, R2initial: Initial concentrations of resources
     Y1, Y2: Yields of resources

     return: Total biomass yield
     """

     # Calculate Rstar
     R1star, R2star = CalcRStarAdditiveModel(gmax, K1, K2, d, R1initial, R2initial, Y1, Y2)

     # Calculate total yield
     dN = (R1initial - R1star)*Y1

     return dN


def CalcResourceDepletionPhaseLine(dR2dR1, params, R1max):
     """
     dR2dR1: Slope of ODE for resource depletion.  Arguments must be R1, R2, and
           *params
     R1max: Maximum value of R1 to end phase line

     return: Interpolated function (R2 as a function of R1) separating resource
          depletion phases
     """

     # Number of points along phase line for ODE calculation, to be used as
     # basis for interpolation
     num_points = 100

     # Minimum values of R1 and R2 to start integration of depletion trajectory
     # --- should be close to zero
     R1min, R2min = 1e-6, 1e-6

     # Set of R1 points along phase line
     R1_trajectory = numpy.linspace(R1min, R1max, num_points)

     # Integrate ODE
     solution = integrate.solve_ivp(dR2dR1, (R1min, R1max), [R2min], t_eval=R1_trajectory, args=params)
     R2_trajectory = solution.y[0]

     # Interpolate trajectory into a smooth function
     phase_line = interpolate.interp1d(R1_trajectory, R2_trajectory, kind="cubic")

     return phase_line


def CalcYieldVariableStoichiometry(R1, R2, Y1_function, Y2_function, phase_line, params):
     """
     R1: Concentration of resource 1
     R2: Concentration of resource 2
     Y1_function: Function of intrinsic yield Y1, call signature should be Y1_function(R1, R2, other parameters)
     Y2_function: Function of intrinsic yield Y2, call signature should be Y2_function(R1, R2, other parameters)
     phase_line: R2 as a function of R1 that separates Omega_1 from Omega_2
     params: Other parameters in model

     return: Total biomass yield at the end of growth
     """

     # Number of points along phase line for ODE calculation, to be used as
     # basis for interpolation
     num_points = 100

     # Minimum values of R1 and R2 to start integration of depletion trajectory
     # --- should be close to zero
     R1min, R2min = 1e-6, 1e-6

     # Determine if this point is in Omega_1
     if R2 > phase_line(R1):

          # Define stoichiometry function for ODE
          dR2dR1 = lambda R1, R2, *params: Y1_function(R1, R2, *params)/Y2_function(R1, R2, *params)

          # Integrate ODE along R1 axis
          R1_trajectory = numpy.linspace(R1, R1min, num_points)
          solution = integrate.solve_ivp(dR2dR1, (R1, R1min), [R2], t_eval=R1_trajectory, args=params)
          R2_trajectory = solution.y[0]
          R2_trajectory_interpolated = interpolate.interp1d(R1_trajectory, R2_trajectory, kind="cubic")
          Y1_function_wrapper = lambda R1: Y1_function(R1, R2_trajectory_interpolated(R1), *params)

          # Now integrate intrinsic yield along the resource depletion trajectory to calculate total yield
          dN = integrate.quad(Y1_function_wrapper, R1min, R1)[0]

     # Otherwise, this point is in Omega_2
     else:

          # Define stoichiometry function for ODE
          dR1dR2 = lambda R2, R1, *params: Y2_function(R1, R2, *params)/Y1_function(R1, R2, *params)

          # Integrate ODE along R2 axis
          R2_trajectory = numpy.linspace(R2, R2min, num_points)
          solution = integrate.solve_ivp(dR1dR2, (R2, R2min), [R1], t_eval=R2_trajectory, args=params)
          R1_trajectory = solution.y[0]
          R1_trajectory_interpolated = interpolate.interp1d(R2_trajectory, R1_trajectory, kind="cubic")
          Y2_function_wrapper = lambda R2: Y2_function(R1_trajectory_interpolated(R2), R2, *params)

          # Now integrate intrinsic yield along the resource depletion trajectory to calculate total yield
          dN = integrate.quad(Y2_function_wrapper, R2min, R2)[0]

     return dN