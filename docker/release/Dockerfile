FROM heltai/deal2lkit:release

MAINTAINER luca.heltai@gmail.com

# pi-BEM master image
RUN git clone https://github.com/mathLab/pi-BEM/ &&\
    mkdir pi-BEM/build && cd pi-BEM/build &&\
    cmake -DCMAKE_BUILD_TYPE=Release \
	  -GNinja \
          ../ && \
    ninja -j4 

# Need this for Travis!
USER root
