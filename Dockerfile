# processor-diffcyt: diffcyt differential abundance / differential state tests on
# FlowSOM metacluster tables. R does the work (R/main.R); a small Python shim
# (processor/) runs it once as a container or serves it as a Lambda function.
# Build with --platform linux/amd64 (the Makefile and docker-compose.yml do).
FROM rocker/r-ver:4.6.1@sha256:bafb3b9a72791eeeb810f0bc0472e74c7f4939b20b46c96cb7ebc6838b6a25ef

# Binary packages from Posit Public Package Manager, frozen at one snapshot date:
# CRAN for R 4.6 and Bioconductor 3.23 (the release that goes with R 4.6).
ARG SNAPSHOT=2026-10-01
ARG BIOC_VERSION=3.23
ENV DEBIAN_FRONTEND=noninteractive

# Runtime libraries the binary packages link against (libxml2 for XML, GLPK/GMP for
# igraph, Cairo/PNG/fonts for the png() device), and Python for the Lambda runtime
# interface client. From the pinned Ubuntu 24.04 base's archive.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      libxml2 libglpk40 libgmp10 libpng16-16t64 libcairo2 libfontconfig1 libfreetype6 \
      python3 python3-venv \
 && rm -rf /var/lib/apt/lists/*

RUN echo "options(repos = c(CRAN = 'https://p3m.dev/cran/__linux__/noble/${SNAPSHOT}'), \
                 BioC_mirror = 'https://p3m.dev/bioconductor/__linux__/noble/${SNAPSHOT}', \
                 HTTPUserAgent = sprintf('R/%s R (%s)', getRversion(), paste(getRversion(), R.version['platform'], R.version['arch'], R.version['os'])))" \
      >> /usr/local/lib/R/etc/Rprofile.site \
 && Rscript -e "install.packages(c('BiocManager', 'jsonlite'))" \
 && Rscript -e "BiocManager::install(version = '${BIOC_VERSION}', ask = FALSE, update = FALSE)" \
 && Rscript -e "BiocManager::install('diffcyt', ask = FALSE, update = FALSE)" \
 && Rscript -e "for (p in c('diffcyt', 'edgeR', 'limma', 'lme4', 'jsonlite')) suppressPackageStartupMessages(library(p, character.only = TRUE)); cat('ok\n')" \
 && rm -rf /tmp/* /usr/local/lib/R/site-library/*/help /usr/local/lib/R/site-library/*/doc /usr/local/lib/R/site-library/*/html

# Lambda runtime interface client, in its own venv (Ubuntu's Python is externally managed).
RUN python3 -m venv /opt/venv && /opt/venv/bin/pip install --no-cache-dir awslambdaric==3.1.1

WORKDIR /app
COPY R/ /app/R/
COPY processor/ /app/processor/
COPY tests/ /app/tests/
COPY VERSION /app/VERSION

ENV PATH="/opt/venv/bin:${PATH}" \
    PYTHONPATH="/app" \
    PYTHONUNBUFFERED=1 \
    APP_ROOT=/app \
    HOME=/tmp

# Runs as uid 1000 on the platform; nothing here needs root.
USER 1000

# Runs R/main.R once as a container, or serves processor.handler as a Lambda
# function (see processor/entry.py).
ENTRYPOINT ["python3", "-m", "processor.entry"]
