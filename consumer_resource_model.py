import numpy
from scipy import integrate
from scipy import interpolate
from scipy import optimize
from matplotlib import pyplot
import colimitation_models


class ConsumerResourceModel:
     """
     This class defines a consumer-resource model and associated functions for 
     calculating its properties and dynamics.
     """

     def __init__(self, gmaxes, resource_thresholds, production_consumption, growth_rate_models, strain_names=None, resource_names=None):
          """
          gmax: List of maximum growth rates for all strains.  Can be list of 
               scalars (one for each strain), or a matrix (rows are strains, 
               columns are resources)
          resource_thresholds: Matrix of resource thresholds (rows are strains,
               columns are resources)
          production_consumption: Matrix of net resource production/consumption 
               per unit biomass (rows are stains, columns are resources).  
               Positive values mean a net production, negative values mean a net 
               consumption.
          growth_rate_models: List of models (one for each strain) to use for
               calculating growth rate with multiple resources
          strain_names: List of strain names
          resource_names: List of resource names
          """

          # Set internal parameters
          self.gmaxes = gmaxes
          self.resource_thresholds = resource_thresholds
          self.production_consumption = production_consumption
          self.growth_rate_models = growth_rate_models

          # Number of strains and resources
          self.num_strains, self.num_resources = self.resource_thresholds.shape

          # If no strain names given, label with alphabet
          if strain_names == None:
               self.strain_names = [chr(65 + i) for i in range(self.num_strains)]
          else:
               if len(strain_names) != self.num_strains:
                    raise ValueError("Inconsistent number of strain names")
               else:
                    self.strain_names = strain_names

          # If no resource names given, label with numbers
          if resource_names == None:
               self.resource_names = [str(i + 1) for i in range(self.num_resources)]
          else:
               if len(resource_names) != self.resource_names:
                    raise ValueError("Inconsistent number of resource names")
               else:
                    self.resource_names = resource_names


     def CalcGrowthRates(self, Rs):
          """
          Rs: List of resource concentrations

          return: Instantaneous growth rates for all strains
          """

          return numpy.array([colimitation_models.CalcTraitODE(Rs, self.gmaxes[s]/self.resource_thresholds[s], self.gmaxes[s], self.growth_rate_models[s]) for s in range(self.num_strains)])


     def CalcDerivatives(self, time, concentrations):
          """
          time: Time value (this is a dummy variable that is not used because 
               the ODE is autonomous, but it is necessary to use the ODE
               integration function)
          concentrations: List of strain and resource concentrations concatenated

          return: Instantaneous derivatives of strain and resource 
               concentrations, concatenated
          """

          # Separate concentrations of strains and resources
          Ns = concentrations[:self.num_strains]
          Rs = concentrations[self.num_strains:]

          # Calculate instantaneous growth rates
          gs = self.CalcGrowthRates(Rs)

          # Calculate strain concentration derivatives
          dNdts = gs*Ns

          # Calculate resource concentration derivatives
          dRdts = dNdts.dot(self.production_consumption)

          return numpy.concatenate([dNdts, dRdts])


     def CalcBatchDynamics(self, N0s, R0s, time_points):
          """
          N0s: Initial concentrations of all strains
          R0s: Initial concentrations of all resources
          time_points: Times at which to calculate concentrations

          return: Separate arrays of strain and resource concentrations at all
               time points
          """
          
          # Integrate ODEs at designated time points
          solution = integrate.solve_ivp(self.CalcDerivatives, (min(time_points), max(time_points)), numpy.concatenate([N0s, R0s]), t_eval=time_points)

          # Separate strain and resource concentrations
          Ns = solution.y[:self.num_strains]
          Rs = solution.y[self.num_strains:]
          return Ns, Rs


     def CalcSaturationTime(self, N0s, R0s, precision=1e-3, mode="convergence"):
          """
          N0s: Initial concentrations of all strains
          R0s: Initial concentrations of all resources
          precision: Relative precision to test strain saturation
          mode: Either "convergence" (standard) or "fold-change" (only valid
               when all nutrients are externally supplied and yields are equal
               across strains)

          return: Time at which strain concentrations appear to have saturated
          """

          if mode == "convergence":                    

               # Relative increase for time step 
               tsat_fold_step = 0.1

               # Initial guess for saturation time; ideally make this more sophisticated so it doesn't overshoot initially
               tsat = 1

               # Initial guess for strain concentrations is their initial values
               Ns_cur = N0s
               Rs_cur = R0s

               # Continue iterating until saturation is reached according to precision parameter
               while True:

                    # Calculate ODE solution up to current tsat guess
                    solution = integrate.solve_ivp(self.CalcDerivatives, (0, tsat), numpy.concatenate([N0s, R0s]), t_eval=(0, tsat))

                    # Get strain concentrations at current tsat guess
                    Ns_next = solution.y[:self.num_strains].T[-1]

                    # If the relative change in strain concentrations is less than the precision, return the current time as saturation time
                    if sum(abs(Ns_next[i] - Ns_cur[i])/Ns_cur[i] for i in range(len(Ns_cur))) < precision:
                         return tsat
                    # Otherwise, increase the time and integrate further
                    else:
                         tsat *= (1 + tsat_fold_step)
                         Ns_cur = numpy.array([n for n in Ns_next])

          elif mode == "fold-change":

               # Get yields for all nutrients --- only for the first strain, and assuming only external nutrients (no secretion)!
               Ys = -1/self.production_consumption[0]

               # Calculate total final biomass, based on the yields
               Nmax = min(R0s*Ys) + sum(N0s)

               # Make an initial guess for the saturation time, if the population grows at the initial rate over the whole fold-change
               tsat = numpy.log(Nmax/sum(N0s))/min(self.CalcGrowthRates(R0s))

               # Now test saturation at the predicted fold-change, and increase time accordingly
               while True:
                    # Calculate ODE solution up to current tsat guess
                    solution = integrate.solve_ivp(self.CalcDerivatives, (0, tsat), numpy.concatenate([N0s, R0s]), t_eval=(0, tsat))

                    # Get strain concentrations at current tsat guess
                    Ns_next = solution.y[:self.num_strains].T[-1]

                    # Calculate total biomass at this time point
                    Ntotal = sum(Ns_next)

                    # If the total biomass at this point is within tolerance of the predicted total, then return the time as tsat
                    if (abs(Ntotal - Nmax)/Nmax < precision) or (Ntotal > Nmax):
                         return tsat
                    tsat *= 2

          else:
               print("Invalid mode for saturation time calculation.  Exiting.")
               exit()


     def CalcSaturationConcentrations(self, N0s, R0s, precision=1e-3, mode="convergence"):
          """
          N0s: Initial concentrations of all strains
          R0s: Initial concentrations of all resources
          precision: Relative precision to test strain saturation

          return: Concentrations of strains at which they appear to have 
               saturated
          """

          # Calculate saturation time
          tsat = self.CalcSaturationTime(N0s, R0s, precision, mode)

          # Solve ODE up to tsat
          solution = integrate.solve_ivp(self.CalcDerivatives, (0, tsat), numpy.concatenate([N0s, R0s]), t_eval=(0, tsat))

          # Get strain concentrations at tsat
          Ns_sat = solution.y[:self.num_strains].T[-1]
          return Ns_sat


     def CalcSelectionCoefficients(self, N0s, R0s, precision=1e-3, mode="convergence"):
          """
          N0s: Initial concentrations of all strains
          R0s: Initial concentrations of all resources
          precision: Relative precision to test strain saturation

          return: Matrix of selection coefficients between all strains at 
               saturation
          """

          # Calculate strain concentrations at saturation
          Ns_sat = self.CalcSaturationConcentrations(N0s, R0s, precision, mode)

          # Calculate selection coefficient matrix
          smatrix = numpy.array([[numpy.log(Ns_sat[i]/Ns_sat[j]) - numpy.log(N0s[i]/N0s[j]) for j in range(self.num_strains)] for i in range(self.num_strains)])

          return smatrix


     # This plots the batch dynamics
     def PlotBatchDynamics(self, time_points, N_trajectories, R_trajectories, save="", precision=1e-3, mode="convergence"):
          """
          time_points: Time points at which to plot the dynamics
          N_trajectories: List of trajectories of strain concentrations
          R_trajectories: List of trajectories of resource concentrations
          save: File name if saving (default is to show)
          """

          # Calculate saturation time
          tsat = self.CalcSaturationTime(N_trajectories.T[0], R_trajectories.T[0], precision, mode)

          # Initialize figure
          figure = pyplot.figure(figsize=(4*2, 3))
          figure.subplots_adjust(hspace=0.3, wspace=0.3)

          # Plot strain concentration trajectories
          axis = figure.add_subplot(1, 2, 1)
          for i in range(len(N_trajectories)):
               axis.plot(time_points, N_trajectories[i], "-", label="Strain " + self.strain_names[i])
          axis.set_xlabel("Time $t$")
          axis.set_ylabel(r"Strain concentration $N_\alpha(t)$")
          axis.set_yscale("log")
          axis.set_xlim([min(time_points), max(time_points)])
          axis.set_ylim([0.9*min(N_trajectories.flatten()), 1.1*max(N_trajectories.flatten())])
          axis.legend(loc="best")
          axis.axvline(tsat, linestyle="--", color="black")

          # Plot resource concentration trajectories
          axis = figure.add_subplot(1, 2, 2)
          for i in range(len(R_trajectories)):
               axis.plot(time_points, R_trajectories[i], "-", label="Resource " + self.resource_names[i])
          axis.set_xlabel("Time $t$")
          axis.set_ylabel("Resource concentration $R_i(t)$")
          axis.set_yscale("log")
          axis.set_xlim([min(time_points), max(time_points)])
          axis.set_ylim([0.9*min(R_trajectories.flatten()), 1.1*max(R_trajectories.flatten())])
          axis.legend(loc="best")
          axis.axvline(tsat, linestyle="--", color="black")

          # Save figure or show
          if save != "":
               figure.savefig(save, bbox_inches="tight")
          else:
               pyplot.show()


     # Calculate trajectories of biomass and nutrients over serial transfers
     def CalcSerialTransfers(self, N0s, R0s, Rsources, time_points_within_cycle, num_cycles, dilution_factor):

          # Initialize biomass and nutrient concentrations
          N0s_cur = numpy.array([N0 for N0 in N0s])
          R0s_cur = numpy.array([R0 for R0 in R0s])

          # Set of all time points
          time_points_over_cycles = []

          # Trajectories (shape is strain, time points)
          Ns_trajectory = numpy.array([[] for N0 in N0s])
          Rs_trajectory = numpy.array([[] for R0 in R0s])

          # Iterate over cycles
          for c in range(num_cycles):

               # Calculate batch dynamics for this cycle
               Ns, Rs = self.CalcBatchDynamics(N0s_cur, R0s_cur, time_points_within_cycle)

               # Dilute biomass
               N0s_cur = numpy.array([N0/dilution_factor for N0 in Ns.T[-1]])

               # Dilute nutrient concentrations and add new concentrations
               R0s_cur = numpy.array([Rs.T[-1][n]/dilution_factor + Rsources[n] for n in range(len(Rs.T[-1]))])

               # Append concentrations from this cycle to ongoing trajectory
               Ns_trajectory = numpy.hstack((Ns_trajectory, Ns))
               Rs_trajectory = numpy.hstack((Rs_trajectory, Rs))

               # Append time points
               time_points_over_cycles = time_points_over_cycles + list(c*max(time_points_within_cycle) + time_points_within_cycle)

          return time_points_over_cycles, Ns_trajectory, Rs_trajectory


     # Calculates the frequency-dependent selection coefficient function for two strains
     def CalcFreqDep(self, N0, R0s, num_freqs=100, min_freq=1e-3, precision=1e-3, mode="convergence"):
          """
          N0: Total initial biomass
          R0s: Initial resource concentrations
          num_freqs: Number of relative frequencies to sample
          min_freq: Minimum frequency to sample

          return: Interpolated function for selection coefficient of second 
               strain over first strain as a function of second strain's initial
               frequency (only valid for two strains)
          """

          # Check that there are only two strains
          if self.num_strains != 2:
               raise ValueError("Cannot calculate frequency dependence for more than two strains")

          # Mesh of initial frequencies for strain 1 and corresponding selection coefficients of strain 1 over strain 2
          x_mesh = numpy.linspace(min_freq, 1 - min_freq, num_freqs)
          selection_mesh = []

          # Iterate over initial frequencies
          for x in x_mesh:

               # Set the initial frequencies
               N0s = N0*numpy.array([1 - x, x])

               # Calculate selection coefficient
               sij = self.CalcSelectionCoefficients(N0s, R0s, precision, mode)

               # Add selection coefficient of strain 2 over strain 1 to mesh
               selection_mesh.append(sij[1][0])

          # Interpolate selection coefficient as a function of strain 2's initial frequency
          interpolated_function = interpolate.interp1d(x_mesh, selection_mesh, kind="cubic")

          return interpolated_function


     # Calculate the coexistence frequency between two strains
     def CalcCoexistenceFreq(self, N0, R0s, num_freqs=100, min_freq=1e-3, precision=1e-3, mode="convergence"):
          """
          N0: Total initial biomass
          R0s: Initial resource concentrations
          num_freqs: Number of relative frequencies to sample
          min_freq: Minimum frequency to sample

          return: Estimated coexistence frequency for second strain over first
               strain (only valid for two strains)
          """

          # Calculate interpolated selection as a function of frequency
          s_of_x = self.CalcFreqDep(N0, R0s, num_freqs, min_freq, precision, mode)

          # Determine the coexistence frequency as the root of the interpolated function
          solution = optimize.root_scalar(s_of_x, bracket=[min_freq, 1 - min_freq])
          return solution.root


     def CalcChemostatDerivatives(self, time, concentrations, Rsources, dilution_rate):

          # Separate concentrations of strains and resources
          Ns = concentrations[:self.num_strains]
          Rs = concentrations[self.num_strains:]

          # Calculate instantaneous growth rates
          gs = self.CalcGrowthRates(Rs)

          # Calculate strain concentration derivatives
          dNdts_growth = gs*Ns
          dNdts_net = (gs - dilution_rate)*Ns

          # Calculate resource concentration derivatives
          dRdts = dNdts_growth.dot(self.production_consumption) + dilution_rate*(Rsources - Rs)

          return numpy.concatenate([dNdts_net, dRdts])


     def CalcChemostatDynamics(self, N0s, R0s, time_points, Rsources, dilution_rate):

          # Integrate ODEs at designated time points
          solution = integrate.solve_ivp(self.CalcChemostatDerivatives, (min(time_points), max(time_points)), numpy.concatenate([N0s, R0s]), t_eval=time_points, args=(Rsources, dilution_rate))

          # Separate strain and resource concentrations
          Ns = solution.y[:self.num_strains]
          Rs = solution.y[self.num_strains:]
          return Ns, Rs

