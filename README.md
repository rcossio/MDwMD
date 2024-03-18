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
- **Use cases**: Identify concrete use cases of the diffusion coefficient of proteins so they can be presented to the community.

### Database

- **Load Data**: Follow HydroPro text to load enough diffusion data and linked PDBs so that you can make a good model.
- **Duplicated uniprot**: Trypsin and trypsinogen have the same uniprot id, as well as chemotrypsin/chemotrypsinogen. We need a new identifier.

### Back-end

- **Good practices**: The project started with a MVP, but we need to apply good programming practices in the backend to have a modular and scalable back-end.
- **Upload to AWS**: Use AWS EC2
- **Conect References colection**: When posting a new experiment we should fetch the reference data

### Front-end

- **Mobile access**: Implementing mobile access would increase productivity (data could be uploaded from the mobile phone, which is easier to use in the tram)
- **Result Sorting**: Fix the sorting of results by accession number, ensuring they appear in order despite varying fetch times.
- **Loading animation and clean-up**: Upon pressing the Search button, it should display an animation of loading, instead of just keeping the old view
- **Data loading**: When loading a new data, trim the information, otherwise wierd character may be copied

### HYDROPRO

- **PDB preparation**: Since HYDROPRO counts non-H atoms from HETATM, most common co-crystalised molecules should be filtered out
- **Crystal waters**: Crystal waters could be part of the system or not, hence we can leave them of remove them. Both situations must be tested to see wich lead to less error.
- **Molecular weight estimation**: (optional) The script that estimates molecular weight gives different results than PDB, it would be nice to reproduce what they calculated.

### Validation

- **Error estimation**: Right now each complex (with possibly multiple experiments) have a single PDB linked. We need to calculate the error associated to each PDB (hence we must treat different PDBs as different samples for the error)
- **Decide experiment outliers**: Some experiments might have outliers. We must stablish criteria to remove them from the DB.
- **Data quantity**: We may have enough data to have a respectable validation set. How much is enough?

### Molecular dynamics
- **Ligandless complexes**: Start by complex without ligand, they are easier


## Insights from Other DBs

- **Data scattering**: "Data is scattered and in that way is useless", a phase I liked.
- **Context**: Provides contextual information about experiments, including pH levels, concentration and temperature.
- **Updates** We should displaying the latest update status prominently on the main page to inform users of the current data accuracy.
- **Experimental diversity**: Covers a variety of experimental conditions, such as different pH levels and protein families, we can focus on protein families, species, molecular weights and distribution of diffusion coefficients. I dont feel confident to include methods, since i dont understand them.
- **Promoting diffusion**: The paper identifies diffusion as the least reported feature, its and opportunity to explain why our DB is relevant.
- **CSV**: we should add the option to download data in CSV format.
- **Data echoing**: We can include their data to make a whole DB (i.e. the data is included as secondary references, or curated to first). Is this scriptable? We can also add the data from Bionumbers with the indexed terms "diffusion", and "diffusion coefficient".
- **Data collection and curation**: Explain how data was treated to show it is trustable. Explain the main keywords used in google scholar searchs. MPTherm explains this well
- **Diluted diffusion**: We have to clearly state that the DB is about diluted solution diffusion (<10mg/ml), preferably infinite diffusion. 
- **Navbars**: Home (overview, updates, cite us), Browse, DB statistics, Tutorial, Upload data   



## Conclusion

This project represents a foundational step towards simplifying the acquisition and estimation of self-diffusion coefficients for proteins, offering significant potential for scientific research and analysis. Through ongoing development and community collaboration, we aim to enhance its utility and impact within the scientific community.
