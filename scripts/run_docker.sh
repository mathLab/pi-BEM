#!/bin/bash
IMG=`grep image .travis.yml | awk '{print $2}'` 
echo $IMG

CMD=`grep -e '- ' .travis.yml | sed 's/- //'`
docker pull $IMG 
docker run -v `pwd`/..:/builds/studenti/pi-BEM/ $IMG /bin/sh -c "$CMD" 
