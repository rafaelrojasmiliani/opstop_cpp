ARG BASEIMAGE
FROM ${BASEIMAGE}
ARG ROS_DISTRO

RUN --mount=type=bind,source=./,target=/workspace,rw \
    cd /workspace \
    && mkdir build \
    && cd build \
    && source /opt/ros/${ROS_DISTRO}/setup.bash \
    && export PKG_CONFIG_PATH=/opt/openrobots/lib/pkgconfig:$PKG_CONFIG_PATH \
    && export LD_LIBRARY_PATH=/opt/openrobots/lib:$LD_LIBRARY_PATH \
    && export CMAKE_PREFIX_PATH=/opt/openrobots:$CMAKE_PREFIX_PATH \
    && cmake .. -DCMAKE_INSTALL_PREFIX=/usr  -DBUILD_TESTING=OFF \
    && make opstop -j $(nproc) \
    && make install
