# Ubuntu 24.04 base image
FROM ubuntu:24.04

# Prevent interactive prompts during package install
ENV DEBIAN_FRONTEND=noninteractive

# Update system and install scientific dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    build-essential \
    gcc \
    g++ \
    make \
    cmake \
    pkg-config \
    libgsl-dev \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Mark mounted repo as safe
RUN git config --global --add safe.directory /workdir

# Create application directory
WORKDIR /workdir

# Copy requirements first for Docker layer caching
COPY requirements.txt .

# Install Python libraries
RUN pip3 install --break-system-packages --no-cache-dir -r requirements.txt

# Default command
CMD ["bash"]