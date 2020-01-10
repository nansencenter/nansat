#!/bin/bash

docker login -u "${DOCKER_USER}" -p "${DOCKER_PASS}"

for image in $*;do
    echo "Tag ${DOCKER_ORG}/${image}:${DOCKER_TAG} as ${DOCKER_ORG}/${image}:latest"
    docker tag "${DOCKER_ORG}/${image}:${DOCKER_TAG}" "${DOCKER_ORG}/${image}:latest"
    echo 'Push the new tags'
    docker push "${DOCKER_ORG}/${image}:${DOCKER_TAG}"
    docker push "${DOCKER_ORG}/${image}:latest"
done