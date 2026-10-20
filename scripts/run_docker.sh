#!/bin/sh
# Run the same warning-free Debug/Release checks used for releases.
exec docker run --rm -v "$PWD:/app" -w /app dealii/dealii:v9.7.1-noble \
  bash -lc 'set -e; cmake -S . -B build-docker/9.7.1 -GNinja -DCMAKE_BUILD_TYPE=DebugRelease -DCMAKE_CXX_FLAGS=-Werror; cmake --build build-docker/9.7.1 -j 4; cd build-docker/9.7.1; ctest --output-on-failure'
