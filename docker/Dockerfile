# docker/Dockerfile
# Dockerfile for Full Spectrum Biomedical Data Architecture

FROM rocker/r-ver:4.3.0

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    git \
    wget \
    unzip \
    curl \
    gsl-bin \
    libgsl-dev \
    python3 \
    python3-pip \
    python3-setuptools \
    bcftools \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# Install PLINK 1.9
RUN mkdir -p /opt/plink && \
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip -O /opt/plink/plink.zip && \
    unzip /opt/plink/plink.zip -d /opt/plink && \
    ln -s /opt/plink/plink /usr/local/bin/plink && \
    rm /opt/plink/plink.zip

# Install Snakemake via pip
RUN pip3 install --break-system-packages snakemake

# Install R Packages
RUN R -e "install.packages(c('plumber', 'RSQLite', 'DBI', 'dplyr', 'ggplot2', 'rmarkdown', 'tinytex'), repos='https://cloud.r-project.org')"

# Setup workspace
WORKDIR /app

# Copy project files
COPY . /app

# Create folder for bin (containing local plink for host execution compatibility)
RUN mkdir -p bin && cp /usr/local/bin/plink bin/plink

# Expose Plumber API port
EXPOSE 8085

# Set default entrypoint to run the API
CMD ["Rscript", "api/app.R"]
