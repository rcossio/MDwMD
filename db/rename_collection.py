from pymongo import MongoClient
import os
from dotenv import load_dotenv

# Load environment variables from a .env file
load_dotenv()
mongo_uri = os.getenv('MONGO_URI')

# Database and collection names
db_name = 'diffusion_db'
original_collection_name = 'complexes'
backup_collection_name = original_collection_name + '_backup'
new_collection_name = 'proteins'

# Connect to MongoDB
client = MongoClient(mongo_uri)
db = client[db_name]

# Check if backup collection already exists to avoid overwriting data
if backup_collection_name in db.list_collection_names():
    raise Exception(f"Backup collection '{backup_collection_name}' already exists. Please choose a different name.")

# Copy documents from the original collection to the backup collection
original_collection = db[original_collection_name]
backup_collection = db[backup_collection_name]

for document in original_collection.find():
    backup_collection.insert_one(document)

print(f"Backup of '{original_collection_name}' completed as '{backup_collection_name}'.")

# Rename the original collection
db[original_collection_name].rename(new_collection_name)

print(f"Collection '{original_collection_name}' has been renamed to '{new_collection_name}'.")
