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
- **hydropro**: Houses scripts for running HYDROPRO calculations.
- **md**: Contains scripts for conducting molecular dynamics simulations.

## Next Steps 

### Concepts
- **Use cases**: Identity use cases of the diffucion coefficient of proteins so they can be presented to the community

### Database

- **PMID Information Table**: Implement a static collection containing PMID information to bypass slow API request times due to frequency limitations of NCBI e-utils.

### Back-end

- **Good practices**: The project started with a MVP, but we need to apply good programming practices in the backend to have a modular and scalable back-end.
- **Upload to AWS**: Use AWS EC2

### Front-end

- **Data Loading and Filtering**: Address the issue of data loading and filtering where old entries interfere with new, filtered entries. Implementing `AbortController` is proposed as a solution.
- **Result Sorting**: Fix the sorting of results by accession number, ensuring they appear in order despite varying fetch times.
- **Mobile access**: Implementing mobile access would increase productivity (data could be uploaded from the mobile phone, which is easier to use in the tram)
- **PDB linking**: There is no current way of linking the PDBs to the complexes with a front-end. Think of a way.

### HYDROPRO

- **PDB preparation**: Since HYDROPRO counts non-H atoms from HETATM, most common co-crystalised molecules should be filtered out
- **Crystal waters**: Crystal waters could be part of the system or not, hence we can leave them of remove them. Both situations must be tested to see wich lead to less error.
- **Molecular weight estimation**: (optional) The script that estimates molecular weight gives different results than PDB, it would be nice to reproduce what they calculated.

### Validation

- **Error estimation**: Right now each complex (with possibly multiple experiments) have a single PDB linked. We need to calculate the error associated to each PDB (hence we must treat different PDBs as different samples for the error)
- **Decide experiment outliers**: Some experiments might have outliers. We must stablish criteria to remove them from the DB.
- **Data quantity**: We don't have enough data to have a respectable validation set.

## Conclusion

This project represents a foundational step towards simplifying the acquisition and estimation of self-diffusion coefficients for proteins, offering significant potential for scientific research and analysis. Through ongoing development and community collaboration, we aim to enhance its utility and impact within the scientific community.
