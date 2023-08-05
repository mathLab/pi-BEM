#!/bin/bash
docker run  --user $(id -u):$(id -g) -i -t --rm -P -v `pwd`:/app:rw dealii/dealii:v9.5.0-jammy  /bin/sh -c "cd /app; $@"
