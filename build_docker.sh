V=`cat ./VERSION`
docker build -f docker/Dockerfile . -t elasticqtl:latest -t elasticqtl:${V}
