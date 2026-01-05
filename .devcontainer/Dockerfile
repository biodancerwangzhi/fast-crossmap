# FastCrossMap Development Container
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install basic system tools
RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    cmake \
    git \
    curl \
    wget \
    vim \
    ca-certificates \
    valgrind \
    heaptrack \
    libclang-dev \
    llvm-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy pre-downloaded Miniforge installer
COPY miniforge.sh /tmp/miniforge.sh

# Install Miniforge
RUN bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

ENV PATH=/opt/conda/bin:$PATH

# Configure conda to use Tsinghua mirrors
RUN conda config --add channels https://mirror.nju.edu.cn/anaconda/pkgs/main && \
    conda config --add channels https://mirror.nju.edu.cn/anaconda/cloud/conda-forge && \
    conda config --add channels https://mirror.nju.edu.cn/anaconda/cloud/bioconda && \
    conda config --set show_channel_urls yes

# Install bioinformatics tools
RUN conda install -y \
    rust \
    htslib \
    samtools \
    python=3.11 \
    crossmap \
    pysam \
    pybigwig \
    bx-python \
    hyperfine \
    fastremap-bio \
    ucsc-liftover \
    && conda clean -afy

# Install Python packages for paper benchmarking and plotting
RUN pip install --no-cache-dir \
    matplotlib \
    numpy \
    pandas \
    seaborn \
    scipy \
    psutil

WORKDIR /workspace

RUN /opt/conda/bin/conda init bash
