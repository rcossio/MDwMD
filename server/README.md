# PDDB web application

## Overview
This repository contains the backend code for the PDDB web application, written in Node.js using the Express framework. The backend is designed to serve HTML files stored in the /public directory. These files contain all the necessary logic as there is no separate frontend build or deployment process.

## Environment Setup
The application uses different environment variables for development and production environments. These variables are specified in the ´.env.development´, ´.env.production´, and ´.env.mongo´ files for the Docker MongoDB container. If you do not have these files, please request them. The application reads these environment variables using the dotenv library.

The default database is hosted on MongoAtlas. For local deployments using Docker, docker-compose sets up a fresh MongoDB instance which is initially empty. Future updates will include configurations to connect automatically to the MongoDB instance created by Docker Compose.

## Usage
Build image and run docker-compose:
```bash
docker compose -p pddb up -d
```

Remove (and delete volumes):
```bash
docker compose -p pddb down #--volumes
```

Enter DB as admin
```bash
mongosh "mongodb://<user>:<pass>@localhost:27017/pddb?authSource=admin"
```
