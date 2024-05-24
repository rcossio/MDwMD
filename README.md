# Project Overview

This project aims to develop a database designed to simplify the process of obtaining and estimating self-diffusion coefficients, which includes validating a prediction model, for a wide range of proteins. By leveraging a comprehensive collection of experimental molecular diffusion data, integrating additional information through APIs like Uniprot and NCBI e-utils, and predicted molecular diffusion values generated through HydroPRO and molecular dynamics simulations we aim to provide a robust resource for researchers and scientists.

## Database Structure

The database contains the following collections:

- **Proteins**: It relates the protein at a certain state or condition with its related PDB structures and experiments. 

- **Experimental Molecular Diffusion**: Contains experimentally measured diffusion coefficients,it can be enhanced when querying with additional data sourced from Uniprot. 

- **References**: It contains reference identifiers, such as PubmedID or DOI. It can be automatically re deployed be means of NCBI e-utils to gather more information from the references. Right now it just contains the year.

- **HydroPRO Predicted Molecular Diffusion**: Includes predictions from HydroPRO, used to validate the method and calculate the error, allowing data extrapolation.


## Project Organization

The project is organized into the following directories, each containing relevant files and scripts. Between brackets who can understand the files within or the main language of the folder:

- **server**: Backend of the web application (Damiano C, Ivan)
- **db**: Constains scripts to mantain and check the DB (python)
- **hydropro**: Has scripts for running HYDROPRO calculations (bash)
- **md**: Contains scripts for conducting molecular dynamics simulations (md: Gonzalo and Giorgia, diffusivity:Franco)
- **validation**: Scripts to evaluate the error of the modelling strategies (python)


## Next Steps 

### Database
- **Load Data**: Follow HydroPro text to load enough diffusion data and linked PDBs so that you can make a good model.
- **MIADE and Disprot**: Read about MIADE and (10.1101/2022.07.12.495092v1) and best practices (https://arxiv.org/pdf/2310.16716.pdf)
- **Alphafold models**: Enable alphafold models as an array attribute in Protein collection

### Back-end
- **Good practices**: The project started with a MVP, but we need to apply good programming practices in the backend to have a modular and scalable back-end.
- **Conect References colection**: When posting a new experiment we should fetch the reference data
- **Filtering negative PDBs**: When displaying Unkown PDBs, we should filter out those that were filtered and discarded

### Front-end
- **Mobile access**: Implementing mobile access would increase productivity (data could be uploaded from the mobile phone, which is easier to use in the tram)
- **Data loading**: When loading a new data, trim the information, otherwise a wierd character may be copied
- **Alphafold**: display alphafold models with the Uniprot API when browsing

### HYDROPRO
- **PDB preparation**: Since HYDROPRO counts non-H atoms from HETATM, most common co-crystalised molecules should be filtered out
- **Crystal waters**: Crystal waters could be part of the system or not, hence we can leave them of remove them. Both situations must be tested to see wich lead to less error.
- **Molecular weight estimation**: (optional) The script that estimates molecular weight gives different results than PDB, it would be nice to reproduce what they calculated.

### Validation
- **Decide experiment outliers**: Some experiments might have outliers. We must stablish criteria to remove them from the DB.
- **Data quantity**: We may have enough data to have a respectable validation set. How much is enough?

### Models
- **Crysol**: Think of using spherical harmonic coefficients to calculate the diffusion coefficient, instead of bead/shell modelling

### Molecular dynamics
- **Diffusivity consistency**: Investigate if the diffusivity valeus achieved are the same for two MDs of the same PDB, and two PDBs of the same protein.


## Insights from Other DBs

- **Data scattering**: "Data is scattered and in that way is useless", a santence I liked.
- **Context**: Provides contextual information about experiments, including pH levels, concentration and temperature.
- **Updates** We should displaying the latest update status prominently on the main page to inform users of the current data accuracy.
- **Experimental diversity**: Covers a variety of experimental conditions, such as different pH levels and protein families, we can focus on protein families, species, molecular weights and distribution of diffusion coefficients. I dont feel confident to include methods, since i dont understand them.
- **Promoting diffusion**: The paper of BioThermDB identifies diffusion as the least reported feature, its and opportunity to explain why our DB is relevant.
- **CSV**: we should add the option to download data in CSV format.
- **Data echoing**: We can include their data to make a whole DB (i.e. the data is included as secondary references, or curated to first). Is this scriptable? We can also add the data from Bionumbers with the indexed terms "diffusion", and "diffusion coefficient".
- **Data collection and curation**: Explain how data was treated to show it is trustable. Explain the main keywords used in google scholar searchs. MPTherm explains this well
- **Diluted diffusion**: We have to clearly state that the DB is about diluted solution diffusion (<10mg/ml), preferably infinite diffusion. 
- **Navbars**: Home (overview, updates, cite us), Browse, DB statistics, Tutorial, Upload data   



## Conclusion

This project represents a foundational step towards simplifying the acquisition and estimation of self-diffusion coefficients for proteins, offering significant potential for scientific research and analysis. Through ongoing development and community collaboration, we aim to enhance its utility and impact within the scientific community.
For any questions please open an issue or write to rodrigoperez93@gmail.com.
