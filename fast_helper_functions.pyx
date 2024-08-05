"""
This module contains two functions to accelerate calculations of the PAT growth rate model.
Compile using python3 setup_fast_helper_functions.py build_ext --inplace
"""

import numpy
from itertools import combinations


def CalcPATSum(double[:] rs):
     """
     rs: List of resource concentrations R divided by their Monod thresholds K

     returns: Sum of terms that goes into the denominator of PAT growth rate
     """

     # Declare total number of resources
     cdef int num_resources = len(rs)

     # Declare and initialize total sum that will be calculated and returned
     cdef double total = 0

     # Declare variable for size of resource subsets
     cdef int resource_subset_size

     # Declare iterator variable
     cdef int i

     # Declare temporary sum variable
     cdef double temp

     # Iterate over number of resources in a subset
     for resource_subset_size in range(1, num_resources + 1):

          # Iterate over combinations of resources for that subset
          for combo in combinations(range(num_resources), resource_subset_size):

               # Initialize temporary sum variable
               temp = 0

               # Sum over resources in the subset
               for i in combo:
                    temp += rs[i]

               # Add contribution from that subset to total sum
               total += (-1)**(resource_subset_size + 1)/temp

     return total


def CalcRateLimitationPAT(int resource_index, double[:] rs):
     # Declare total number of resources
     cdef int num_resources = len(rs)

     # Growth rate divided by gmax
     cdef double growth_rate_divided_gmax = 1/(1 + CalcPATSum(rs))

     # Declare and initialize total sum that will be calculated and returned
     cdef double Li = 1/rs[resource_index]**2

     # Declare variable for size of resource subsets
     cdef int resource_subset_size

     # Declare iterator variable
     cdef int i

     # Declare temporary sum variable
     cdef double temp

     cdef list indices = [j for j in range(num_resources) if j != resource_index]

     # Iterate over number of resources in a subset
     for resource_subset_size in range(1, num_resources):

          # Iterate over combinations of resources for that subset
          for combo in combinations(indices, resource_subset_size):

               # Initialize temporary sum variable
               temp = rs[resource_index]

               # Sum over resources in the subset
               for i in combo:
                    temp += rs[i]

               # Add contribution from that subset to total sum
               Li += (-1)**resource_subset_size/temp**2

     return Li*rs[resource_index]*growth_rate_divided_gmax