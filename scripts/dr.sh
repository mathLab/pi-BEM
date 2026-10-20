#!/bin/sh
# Run a shell command in the supported deal.II environment.
exec docker run --user "$(id -u):$(id -g)" -it --rm \
  -v "$PWD:/app" -w /app dealii/dealii:v9.8.0-resolute-r1-arm64 /bin/sh -c "$*"
