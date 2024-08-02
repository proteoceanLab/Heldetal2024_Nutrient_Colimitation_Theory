# Heldetal2024_Nutrient_Colimitation_Theory
Code for data analysis and plotting for Held et al., 2024 "Nutrient colimitation is a quantitative, dynamic property of microbial populations"

# Contains:
- /Data: folder containing data files required for plotting
- colimitation_plots.py: general functions for plotting colimitation
- setup_fast_helper_functions.py: helper functions for plotting
- Python scripts for Figures 1-5
- Python script for Figures S1-S31

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

