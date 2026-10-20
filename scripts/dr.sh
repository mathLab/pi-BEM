#!/bin/sh
# Run a shell command in the supported deal.II environment.
exec docker run --user "$(id -u):$(id -g)" -it --rm \
  -v "$PWD:/app" -w /app dealii/dealii:v9.7.1-noble /bin/sh -c "$*"
