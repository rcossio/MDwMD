from pymongo import MongoClient
import json
import os
from dotenv import load_dotenv

def get_db_connection():
    """Establish a connection to the MongoDB database."""
    load_dotenv()  # Load environment variables from a .env file
    uri = os.getenv('MONGO_URI')
    if uri is None:
        raise ValueError("MONGO_URI environment variable not set.")
    client = MongoClient(uri)
    return client['diffusion_db']

def fetch_data(db):
    """Fetch data from MongoDB collections and process it."""
    data = []
    experiments = db['data_from_experiments'].find({
        "active": True,
        "$and": [{"pdbStructures": {"$ne": []}}, {"pdbStructures": {"$ne": ""}}]
    })

    for experiment in experiments:
        process_document(experiment, data)
        print('.')

    return data

def process_document(doc, data):
    """Process a single document from the database."""
    accession_number = doc.get('accessionNumber', '')
    experimental_diffusion_coefficient = get_experimental_diffusion_coefficient(doc)

    pdb_list = doc['pdbStructures']

    update_data_list(accession_number, experimental_diffusion_coefficient, pdb_list, data)

def get_experimental_diffusion_coefficient(doc):
    """Convert diffusion coefficient to a unified scale based on units."""
    unit = doc.get('diffusionUnit')
    coefficient = doc.get('diffusionCoefficient', '')
    if coefficient == '':
        raise ValueError("Missing experimental data for document.")

    if unit == "1e-7 cm2/s":
        return coefficient * 10
    elif unit == "um2/s":
        return coefficient
    else:
        raise ValueError("Unknown unit of measurement.")

def update_data_list(accession_number, experimental_diffusion_coefficient, pdb_list, data):
    """Update the data list with the processed information."""
    for item in data:
        if item['accessionNumber'] == accession_number:
            item['experimentalDiffusionCoefficient'].append(experimental_diffusion_coefficient)
            item['pdbList'].extend(pdb_list)
            break
    else:  # This else corresponds to the for loop, executed only if the loop completes normally (no break)
        data.append({
            "accessionNumber": accession_number,
            "experimentalDiffusionCoefficient": [experimental_diffusion_coefficient],
            "pdbList": pdb_list
        })

def main():
    """Main function to orchestrate data fetching and processing."""
    db = get_db_connection()
    data = fetch_data(db)
    
    # Iterate over each item in data to update with calculated diffusion coefficients
    for item in data:
        pdb_list = item['pdbList']  # Get the pdbList for the current item
        unique_pdb_codes = set(pdb_list)  # Ensure uniqueness within this list
        
        calculated_diffusion_coefficients = []
        for pdb_code in unique_pdb_codes:
            hydropro_doc = db['hydropro_calcs'].find_one({"pdbCode": pdb_code})
            if hydropro_doc:
                coefficient = hydropro_doc.get('diffusionCoefficient')
                if coefficient is not None:
                    calculated_diffusion_coefficients.append(coefficient)
        
        # Update the item with its calculated diffusion coefficients
        item['calculatedDiffusionCoefficient'] = calculated_diffusion_coefficients

    # Write data as a json file
    with open('diffusion_data.json', 'w') as f:
        json.dump(data, f, indent=4)

if __name__ == "__main__":
    main()