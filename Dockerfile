# Orchestrator-only runtime image.
#
# This image provides locked Python 3.12, Nextflow, Java, Biopython, pandas
# and Typer. BLAST+, DIAMOND, AUGUSTUS, HMMER, InterProScan, Kallisto, Comet,
# Percolator and reference databases are NOT included. Mount them into /data
# or provide them through a site-specific Nextflow configuration.

FROM ghcr.io/prefix-dev/pixi:0.63.2@sha256:034c9032401f39bbc65d4cb533ac4f834df9d994085558cec75ed800e76a9424

LABEL org.opencontainers.image.title="AnnotaPipeline orchestrator"
LABEL org.opencontainers.image.description="Locked Python and Nextflow runtime for AnnotaPipeline; scientific tools are not bundled."

ENV PIXI_HOME=/opt/pixi \
    NXF_HOME=/work/.nextflow \
    PATH=/opt/annotapipeline/.pixi/envs/runtime/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

WORKDIR /opt/annotapipeline

COPY pixi.toml pixi.lock ./
RUN pixi install --locked --environment runtime \
    && pixi clean cache --yes

COPY main.nf nextflow.config LICENSE ./
COPY lib ./lib
COPY modules ./modules
COPY subworkflows ./subworkflows
COPY workflows ./workflows
COPY scripts ./scripts
COPY config ./config

RUN mkdir -p /data /work \
    && useradd --create-home --uid 10001 annota \
    && chown -R annota:annota /opt/annotapipeline /data /work

VOLUME ["/data", "/work"]
WORKDIR /work
USER annota

CMD ["bash"]
