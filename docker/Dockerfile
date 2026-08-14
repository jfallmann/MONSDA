FROM mambaorg/micromamba:latest

ARG ENVIRONMENT
ARG MONSDA_VERSION
ARG VCS_REF

LABEL org.opencontainers.image.source="https://github.com/jfallmann/MONSDA"
LABEL org.opencontainers.image.revision="${VCS_REF}"
LABEL org.opencontainers.image.version="${MONSDA_VERSION}"
LABEL org.opencontainers.image.title="MONSDA ${ENVIRONMENT} environment"

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/${ENVIRONMENT}.yaml /tmp/environment.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts /opt/MONSDA/scripts

RUN micromamba create --yes --prefix /opt/conda/envs/monsda --file /tmp/environment.yaml \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/envs/monsda/bin:${PATH}"
ENV MONSDA_VERSION="${MONSDA_VERSION}"

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["/bin/bash"]
