from pymongo.mongo_client import MongoClient
import json

# Establish a connection to the MongoDB database
uri = "mongodb+srv://guest:4321guest@cluster0.hgqsmv1.mongodb.net/?retryWrites=true&w=majority&appName=Cluster0"
# Create a new client and connect to the server
client = MongoClient(uri)
db = client['diffusion_db']

# Access the collections
data_from_experiments = db['data_from_experiments']
hydropro_calcs = db['hydropro_calcs']

# Query the "data_from_experiments" collection
query = {
    "species": "9606",
    "active": True,
    "$and": [
        {"pdbStructures": {"$ne": []}},  # Not an empty array
        {"pdbStructures": {"$ne": ""}}   # Not an empty string
    ]
}
documents = data_from_experiments.find(query)

data = []

# Process each document
for doc in documents:
    accession_number = doc.get('accessionNumber', '')
    for pdb_code in doc['pdbStructures']:
        # Find the related document in "hydropro_calcs" based on pdb_code
        hydropro_doc = hydropro_calcs.find_one({"pdbCode": pdb_code})
        if hydropro_doc:
            calculated_diffusion_coefficient = hydropro_doc.get('diffusionCoefficient', 'N/A')
            if doc.get('diffusionUnit') == "1e-7 cm2/s":
                experimental_diffusion_coefficient = doc.get('diffusionCoefficient', '') * 10
            elif doc.get('diffusionUnit') == "um2/s":
                experimental_diffusion_coefficient = doc.get('diffusionCoefficient', '')
            else:
                experimental_diffusion_coefficient = 'N/A'
            # Find the item in data that has the same accession number and add the experimental diffusion coeficient, if there is none add one entry with all the data 
            for item in data:
                if item['accessionNumber'] == accession_number:
                    item['experimentalDiffusionCoefficient'].append(experimental_diffusion_coefficient)
                    break
            else:
                data.append({
                    "accessionNumber": accession_number,
                    "experimentalDiffusionCoefficient": [experimental_diffusion_coefficient],
                    "calculatedDiffusionCoefficient": calculated_diffusion_coefficient
                })

# Write data as a json file
with open('diffusion_data.json', 'w') as f:
    json.dump(data, f, indent=4)    