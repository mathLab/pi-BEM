#!/bin/bash
IMG=`grep image ../.gitlab-ci.yml | awk '{print $2}'` 
echo $IMG

CMD=`grep -e '- ' ../.gitlab-ci.yml | sed 's/- //'`
docker pull $IMG 
docker run -v `pwd`/..:/builds/studenti/pi-BEM/ $IMG /bin/sh -c "$CMD" 
