#!/bin/bash
IMG=mathlab/deal2lkit:latest #`grep image .travis.yml | awk '{print $2}'` 
echo $IMG

CMD=`grep -e '- ' .travis.yml | sed 's/- //'`
echo $CMD
docker pull $IMG 
docker run -v `pwd`/..:/app/ $IMG /bin/sh -c "cd /app && $CMD" 
