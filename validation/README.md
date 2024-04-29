# Validation

To obtain plots to valdidate the calculation of diffusion coefficients you must:
1. Get the data from the DB (hydropro or MD collections, plus the experimental values). This is done with the script ´1_get_data_<type>.py´ and outputs the data to a json file. 
2. (Optional) Obtain the diffusivity, experimental and predicted, listed on top of the protein labels. This is done with the script ´2_range_of_values.py´.
3. Get the experimental vs predicted plot, with the estimation of the error. This si done with the script ´3_value_vs_value.py´.