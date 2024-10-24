# Heldetal2024_Nutrient_Colimitation_Theory
Code for data analysis and plotting for Held et al., 2024 "Nutrient colimitation is a quantitative, dynamic property of microbial populations"

# Contains:
- /Data: folder containing data files required for plotting
- Analyze_rate_data.py and Analyze_yield_data.py: Fit models and calculate supplementations for growth rate and growth yield glucose-ammonium scans
- Bootstrap_rate_data.py and Bootstrap_yield_data.py: Bootstrap data from growth rate and growth yield glucose-ammonium scans and fit models
- Simulate_rate_data.py and Simulate_yield_data.py: Simulate growth rate and growth yield glucose-ammonium scans using models fit to data, fit models to simulated data, and calculate supplementations
- colimitation_models.py: Functions for models of growth rate and growth yield
- colimitation_data_analysis.py: Functions for analyzing growth rate and growth yield glucose-ammonium scans
- colimitation_plots.py: Functions and settings for making plots
- consumer_resource_model.py: Class for integrating dynamics of a consumer-resource model ODE
- fast_helper_functions.pyx and setup_fast_helper_functions.py: Scripts for Cython code to efficiently calculate growth rate and limitation coefficients for PAT model
- Python scripts for Figures 1-5
- Python scripts for Figures S1-S32

# Requirements

This code requires the following packages, some of which may require updated versions to function correctly:

* Python version 3.10.9
* NumPy version 1.24.1
* SciPy version 1.10.0
* Pandas version 1.5.2
* Matplotlib version 3.9.1
* Cython version 0.29.33

Compile the Cython code fast helper functions (for the PAT/SU model):

```
python3 setup_fast_helper_functions.py build_ext --inplace

