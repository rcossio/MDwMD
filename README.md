# Project Overview

This project aims to develop a database designed to simplify the process of obtaining and estimating self-diffusion coefficients for a wide range of proteins. By leveraging a comprehensive collection of experimental molecular diffusion data and integrating additional information through APIs from Uniprot and NCBI e-utils, we aim to provide a robust resource for researchers and scientists. Furthermore, the project includes collections of predicted molecular diffusion values generated through HydroPRO and molecular dynamics simulations, facilitating the analysis of these methods' predictive capabilities on unknown datasets.

## Database Structure

The database is structured into several key components, each serving a specific function within the broader scope of the project:

- **Experimental Molecular Diffusion**: A collection of measured diffusion coefficients, enhanced with additional data sourced from Uniprot and NCBI e-utils.
- **Predicted Molecular Diffusion**: Includes predictions from HydroPRO and molecular dynamics, aimed at evaluating and extrapolating these methods to uncharted data.

## Project Organization

The project is organized into the following directories, each containing relevant files and scripts:

- **server**: Contains the web application for database access.
- **data**: Reserved for initial data and processing scripts, currently unused.
- **validation**: Scripts to evaluate the performance of the modelling strategies
- **hydropro**: Houses scripts for running HydroPRO calculations.
- **md**: Contains scripts for conducting molecular dynamics simulations.

## Next Steps 

### Database

- **PMID Information Table**: Implement a static collection containing PMID information to bypass slow API request times due to frequency limitations of NCBI e-utils.

### Front-end

- **Data Loading and Filtering**: Address the issue of data loading and filtering where old entries interfere with new, filtered entries. Implementing `AbortController` is proposed as a solution.
- **Result Sorting**: Fix the sorting of results by accession number, ensuring they appear in order despite varying fetch times.

### HYDROPRO

- **PDB preparation**: Since HYDROPRO counts non-H atoms from HETATM, most common co-crystalised molecules should be filtered out

- **Crystal waters**: Crystal waters could be part of the system or not, hence we can leave them of remove them. Both situations must be tested to see wich lead to less error.

## Conclusion

This project represents a foundational step towards simplifying the acquisition and estimation of self-diffusion coefficients for proteins, offering significant potential for scientific research and analysis. Through ongoing development and community collaboration, we aim to enhance its utility and impact within the scientific community.
