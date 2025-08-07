FROM rootproject/root:6.32.00-ubuntu22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    nlohmann-json3-dev \
    libboost-all-dev \
    python3 \
    python3-pip && \
    apt-get clean

# Upgrade pip and install Python packages
RUN pip3 install --upgrade pip
RUN pip3 install numpy scipy matplotlib

WORKDIR /code
COPY . /code
RUN rm -rf /code/build

RUN mkdir -p build && cd build && cmake .. && make -j$(nproc) install

CMD ["/bin/bash"]
