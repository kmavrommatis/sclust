# Building Sclust for Container Deployment

## Overview
This guide explains how to build and install Sclust in a container where the binary will be located in `/usr/local/bin` and support files in standard Unix locations.

## Changes Made

### 1. Configurable INSTALLDIR
The makefile now supports setting `INSTALL_PREFIX` to specify where Sclust will look for R scripts:

```bash
# For local development (default)
make

# For container installation
make INSTALL_PREFIX=/usr/local/share/Sclust
```

### 2. Install Target
A new `install` target copies files to standard locations:
- Binary: `/usr/local/bin/Sclust`
- R scripts: `/usr/local/share/Sclust/R/`

### 3. PATH-based Execution
The code now calls `Sclust cluster` directly instead of using hardcoded paths, allowing it to find the binary in PATH.

## Building for Containers

### Dockerfile Example

```dockerfile
FROM ubuntu:20.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    r-base \
    libz-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy source code
WORKDIR /build
COPY . .

# Build with install prefix for container deployment
RUN make clean && \
    make INSTALL_PREFIX=/usr/local/share/Sclust && \
    make install INSTALL_PREFIX=/usr/local/share/Sclust

# Cleanup build artifacts
WORKDIR /
RUN rm -rf /build

# Set working directory for runtime
WORKDIR /data
```

### Manual Build and Install

```bash
# Clean previous builds
make clean

# Build with container paths
make INSTALL_PREFIX=/usr/local/share/Sclust

# Install to system locations (may require sudo)
sudo make install INSTALL_PREFIX=/usr/local/share/Sclust
```

### Verification

After installation, verify the setup:

```bash
# Check binary
which Sclust
# Should output: /usr/local/bin/Sclust

# Check R scripts
ls /usr/local/share/Sclust/R/
# Should show: plot_cn.R  plot_cluster.R

# Test execution
Sclust --help
```

## Important Notes

### bamprocess Data Directory
The `bamprocess` command requires additional reference data (partitions, annotations). In containers, you have two options:

1. **Mount the data directory**: 
   ```bash
   docker run -v /path/to/data:/data my-sclust-image \
     Sclust bamprocess --dir /data ...
   ```

2. **Copy data into the container**:
   ```dockerfile
   COPY partitions /usr/local/share/Sclust/partitions
   COPY annotation /usr/local/share/Sclust/annotation
   ```
   Then modify code to use `INSTALLDIR` for data_dir default.

### R Plotting
The `cn` command will only generate plots if R is installed. The makefile detects R and sets the `RPLOTTING` flag automatically. Without R, the program runs but skips plot generation.

## Rebuilding installdir.h

The `installdir.h` file is auto-generated during build. If you need to manually regenerate it:

```bash
make installdir.h INSTALL_PREFIX=/usr/local/share/Sclust
```

## Troubleshooting

### Error: "cannot open file .../R/plot_cn.R"
This means INSTALLDIR was set incorrectly during compilation. Rebuild with:
```bash
make clean
make INSTALL_PREFIX=/usr/local/share/Sclust
make install INSTALL_PREFIX=/usr/local/share/Sclust
```

### Error: "Sclust: command not found"
The binary is not in PATH. Either:
- Add `/usr/local/bin` to PATH
- Use full path: `/usr/local/bin/Sclust`
- Ensure `make install` completed successfully
