FROM dealii/dealii:v9.8.0-resolute-r1-arm64

# pi-BEM master image
RUN git clone https://github.com/mathLab/pi-BEM/ &&\
    mkdir pi-BEM/build && cd pi-BEM/build &&\
    cmake -DCMAKE_BUILD_TYPE=DebugRelease \
    -GNinja \
    ../ && \
    ninja
