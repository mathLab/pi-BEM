#!/bin/sh
# Open the supported deal.II environment in the current checkout.
exec docker run --rm -it -v "$PWD:/app" -w /app dealii/dealii:v9.8.0-resolute-r1-arm64 bash
