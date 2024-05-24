from pymongo import MongoClient
import os
from dotenv import load_dotenv

# Load environment variables from a .env file
load_dotenv()
mongo_uri = os.getenv('MONGO_URI')

# Database names
original_db_name = 'diffusion_db'
new_db_name = 'pddb_static_copy'

# Connect to MongoDB
client = MongoClient(mongo_uri)

# Access original database
original_db = client[original_db_name]

# Create or access new database
new_db = client[new_db_name]

# Get list of all collections in the original database
collections = original_db.list_collection_names()

# Copy each collection from the original database to the new database
for collection_name in collections:
    original_collection = original_db[collection_name]
    new_collection = new_db[collection_name]
    
    # Copy all documents from the original collection to the new collection
    for document in original_collection.find():
        new_collection.insert_one(document)

print(f"Database '{original_db_name}' has been copied to '{new_db_name}'.")
