FROM rootproject/root:6.32.00-ubuntu22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    software-properties-common && \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get update && apt-get install -y \
    gcc-13 \
    g++-13 \
    cmake \
    nlohmann-json3-dev \
    libboost-all-dev \
    python3 \
    python3-pip \
    libfmt-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Make GCC 13 the default
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 100 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-13 100

# Upgrade pip and install Python packages
RUN pip3 install --upgrade pip
RUN pip3 install numpy scipy matplotlib

WORKDIR /code
COPY . /code
RUN rm -rf /code/build

RUN mkdir -p build && cd build && cmake .. && make -j$(nproc) install

CMD ["/bin/bash"]
