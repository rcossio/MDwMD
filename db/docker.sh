# Line to run mongo db in docker that sets:
# - 2GB of memory
# - 2 CPUs
# - 10MB max size of log files && 3 max number of log files
# - port 27017 is exposed

docker run -d -p 27017:27017  --memory 2g --cpus=2  --log-opt max-size=10m --log-opt max-file=3  --name mongo-container-0 mongo:latest

# Line to run the backend
docker run -d -p 3000:3000  --memory 2g --cpus=2  --log-opt max-size=10m --log-opt max-file=3  --name pddt-container-0 pddb-backend-img

# Line to create a network
docker network create pddt-network-0
