# syntax=docker/dockerfile:1

ARG SEURAT_VERSION=5.4.0
FROM satijalab/seurat:${SEURAT_VERSION}

LABEL authors="Michael Weinberger <michael.weinberger11@gmail.com>"

# install packages
RUN apt-get update && \
	apt-get install -y libcurl4-openssl-dev libssl-dev procps && \
    apt-get clean -y && \
    rm -rf /var/lib/apt/lists/*

# Install R dependencies
ENV R_DEFAULT_INTERNET_TIMEOUT=3600

RUN install2.r --error \ 
    BiocManager \
    cowplot \
    dplyr \
    ggplot2 \
    gprofiler2 \
    hdf5r \
    Matrix \
    patchwork \
    pheatmap \
    RColorBrewer \
    readr \
    remotes \
    SCINA \
    sctransform \
    SoupX \
    stringr \
    tidyr \
    viridis \
    writexl

RUN R -e "BiocManager::install('anndataR')" \
&& R -e "BiocManager::install('biomaRt')" \
&& R -e "BiocManager::install('BSgenome')" \
&& R -e "BiocManager::install('DropletUtils')" \
&& R -e "BiocManager::install('GenomicFeatures')" \
&& R -e "BiocManager::install('rtracklayer')" \
&& R -e "BiocManager::install('txdbmaker')" \
&& R -e "BiocManager::install('XVector')" \
&& R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')" \
&& R -e "BiocManager::install('BSgenome.Mmusculus.UCSC.mm39')" \
&& R -e "BiocManager::install('BSgenome.Drerio.UCSC.danRer11')"

RUN R -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')" \
&& R -e "remotes::install_github('immunogenomics/harmony')" 

# Create a non-privileged user that the app will run under.
# See https://docs.docker.com/go/dockerfile-user-best-practices/
ARG UID=10001
RUN adduser \
    --disabled-password \
    --gecos "" \
    --home "/nonexistent" \
    --shell "/sbin/nologin" \
    --no-create-home \
    --uid "${UID}" \
    appuser

# Switch to the non-privileged user
USER appuser

# by default /bin/bash is executed
#ENTRYPOINT ["/bin/bash"]