#!/bin/bash

for tag in $*;do
    echo "Tag ${IMAGE_NAME}:${DOCKER_TMP_TAG} as ${IMAGE_NAME}:${tag}"
    docker tag "${IMAGE_NAME}:${DOCKER_TMP_TAG}" "${IMAGE_NAME}:${tag}"
    echo "Push ${IMAGE_NAME}:${tag}"
    docker push "${IMAGE_NAME}:${tag}"
done