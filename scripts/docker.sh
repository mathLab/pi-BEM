#!/bin/sh
# Open the supported deal.II environment in the current checkout.
exec docker run --rm -it -v "$PWD:/app" -w /app dealii/dealii:v9.7.1-noble bash
