#!/bin/sh
# Run the same warning-free Debug/Release checks used for releases.
exec docker run --rm -v "$PWD:/workspace" -w /workspace dealii/dealii:v9.8.0-resolute-r1-arm64 \
  bash -lc 'set -e; cmake -S . -B build-docker/9.8.0 -GNinja -DCMAKE_BUILD_TYPE=DebugRelease -DCMAKE_CXX_FLAGS=-Werror; cmake --build build-docker/9.8.0 -j 4; cd build-docker/9.8.0; ctest --output-on-failure'
