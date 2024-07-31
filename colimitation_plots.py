import numpy
import math
import pandas
import sys
#sys.path.append("/home/mmanhart/Tools")
#import crm
#import colimitation


# Set font sizes
panel_label_size = 20
axis_label_size = 10
tick_label_size = 10
legend_label_size = 8
legend_label_size_small = 6


# Set default gmax
gmaxes = numpy.array([1])


# Set default dilution factor
dilution_factor = 100


# Set colors for nutrients
nutrient_colors = ["blue", "orange", "green", "red", "purple", "brown", "pink", "greenyellow", "olive", "cyan"]


def SetrcParams(rcParams):

     # Set fonts
     rcParams["font.family"] = "serif"
     rcParams["text.usetex"] = True
     rcParams["font.size"] = tick_label_size

     return rcParams


def GetMooreYields(num_nutrients, normalize=True, lim_nutrient=1):

     # Read data table from Moore et al. paper
     data = pandas.read_csv("../Moore_etal_2013_NatureGeosci_data/Supplementary_Table_1.csv")

     # List of pairs: elements and relative biomass yields
     pair_list = []

     # Iterate over rows in data table
     for index, row in data.iterrows():

          # Skip elements with missing data for these entries (H and O)
          if (not math.isnan(row["Phytoplankton quota (mol:mol C)"])) and (not math.isnan(row["Mean ocean concentration (umol/kg)"])):

               # Get relative yield (biomass per unit element)
               Y = 1/row["Phytoplankton quota (mol:mol C)"]

               # Get element concentration and convert to mol from umol
               R = row["Mean ocean concentration (umol/kg)"]/1e6

               # Add data to list
               pair_list.append((row["Element"], R*Y))

     # Sort elements and ratios in order of increasing biomass
     pair_list = sorted(pair_list, key=lambda x: x[1])

     # Extract sorted lists of elements and yields separately; exclude Cd
     elements = [pair[0] for pair in pair_list if pair[0] != "Cd"]
     yields = numpy.array([pair[1] for pair in pair_list if pair[0] != "Cd"])

     # Extract lowest yields and normalize the smallest to 1
     yields = yields[:num_nutrients]
     if normalize:
          yields = yields/yields[0]

     # Reorder yields so that a potentially different yield is the lowest
     yields[0], yields[lim_nutrient] = yields[lim_nutrient], yields[0]

     return yields


def ReadSimulationData(run_name):

     # Parameters
     params = {}

     # Read simulation data
     with gzip.open(run_name + ".K-trajectories.gz", "rb") as myFile:

          # Initialize trajectories data array
          trajectories = []

          # Read lines in file
          for line in myFile:

               # Convert binary to ASCII
               line = line.decode("ascii")

               # Skip blank lines
               if line.isspace(): continue

               # Check for parameter lines
               if line.startswith("#"):

                    # Extract growth rate model and yields
                    if line.startswith("# Namespace"):
                         params["growth_rate_model"] = line[line.find("growth_rate_model='") + len("growth_rate_model='"):line.find("', mean_mutation_effect")]

                         position = line.find("yields=")
                         line = line[position + len("yields="):].strip().strip(")").strip("'").split(",")
                         params["Ys"] = numpy.array([float(x) for x in line])

                    # Extract gmax
                    elif line.startswith("# gmaxes"):
                         line = line.strip().strip("# gmaxes").strip("[").strip("]").split()
                         params["gmaxes"] = numpy.array([numpy.mean([float(x) for x in line])])

                    # Extract initial biomass
                    elif line.startswith("# N0s"):
                         line = line.strip().strip("# N0s").strip("[").strip("]").split()
                         params["N0s"] = numpy.array([sum(float(x) for x in line)])
                    
                    # Extract initial nutrient concentrations
                    elif line.startswith("# R0s"):
                         line = line.strip().strip("# R0s").strip("[").strip("]").split()
                         params["R0s"] = numpy.array([float(x) for x in line])
                         params["num_nutrients"] = len(params["R0s"])

                    continue

               # Split the line into the trajectory data
               line = line.split()
               trajectories.append([[float(K) for K in entry.split(",")] for entry in line])

     # Convert trajectory data to numpy array (shape is replicates, time points, nutrients)
     trajectories = numpy.array(trajectories)
     
     # Average over replicates to get average trajectory
     trajectory_average = numpy.average(trajectories, axis=0)

     # Extract evolutionary time points (number of mutation events)
     time_points_evo = numpy.array([t for t in range(len(trajectory_average))])

     return trajectories, trajectory_average, time_points_evo, params


def ReadMeffTrajectories(run_name):

     time_points = []
     Meff_trajectories = []

     with open(run_name + ".Meff-trajectories") as myFile:
          for line in myFile:
               if line.isspace() or line.startswith("#"): 
                    continue

               time_points.append(int(line.split()[0]))
               Meff_trajectories.append([float(x) for x in line.split()[1:]])

     time_points = numpy.array(time_points)
     Meff_trajectories = numpy.array(Meff_trajectories).T

     return time_points, Meff_trajectories

