from pymongo import MongoClient
import json
import os
from dotenv import load_dotenv
from bson.objectid import ObjectId 

def get_db_connection():
    """Establish a connection to the MongoDB database."""
    load_dotenv()  # Load environment variables from a .env file
    uri = os.getenv('MONGO_URI')
    if uri is None:
        raise ValueError("MONGO_URI environment variable not set.")
    client = MongoClient(uri)
    return client['diffusion_db']

def get_experiment_coef(id, db):
    """Convert diffusion coefficient to a unified scale based on units."""
    experiment = db['data_from_experiments'].find_one({"_id": ObjectId(id)}) 

    unit = experiment.get('diffusionUnit')
    coefficient = experiment.get('diffusionCoefficient', '')

    if unit == "1e-7 cm2/s":
        return coefficient * 10
    elif unit == "um2/s":
        return coefficient
    else:
        raise ValueError("Unknown unit of measurement.")

def get_calculated_coef(pdb_code, db):
    """Get the calculated diffusion coefficient from the hydropro_calcs collection."""
    hydropro_doc = db['hydropro_calcs'].find_one({"pdbCode": pdb_code})
    if hydropro_doc:
        coefficient = hydropro_doc.get('diffusionCoefficient')
        if coefficient is not None:
            return coefficient

def main():
    """Main function to orchestrate data fetching and processing."""
    db = get_db_connection()
    data = []
    complexes = db['complexes'].find({})
    print('Processing', complexes.count(), 'complexes...')

    for doc in complexes:

        item = {
            "id": str(doc['_id']),
            "name": doc['name'],
            "experiment_coef": [],
            "calculated_coef": []
        }

        # Get the experimental values
        exp_list = doc['experiments']
        for exp_id in exp_list:
            coef = get_experiment_coef(exp_id, db)
            item['experiment_coef'].append(coef)

        # Get the calculated values 
        pdb_list = doc['pdbStructures']  # Get the pdbList for the current item
        for pdb_code in pdb_list:
            item['calculated_coef'].append(get_calculated_coef(pdb_code, db))

        data.append(item)

    print(len(data))
    # Write data as a json file
    with open('diffusion_data.json', 'w') as f:
        json.dump(data, f, indent=4)

if __name__ == "__main__":
    main()