FROM condaforge/miniforge3:latest

ARG ENVIRONMENT
ARG MONSDA_VERSION
ARG VCS_REF

LABEL org.opencontainers.image.source="https://github.com/jfallmann/MONSDA"
LABEL org.opencontainers.image.revision="${VCS_REF}"
LABEL org.opencontainers.image.version="${MONSDA_VERSION}"
LABEL org.opencontainers.image.title="MONSDA ${ENVIRONMENT} environment"

RUN useradd --create-home --uid 57439 monsda

COPY envs/${ENVIRONMENT}.yaml /tmp/environment.yaml
COPY scripts /opt/MONSDA/scripts
RUN find /opt/MONSDA/scripts -type f \( -name '*.py' -o -name '*.R' -o -name '*.sh' -o -name '*.pl' \) -exec chmod a+x {} +

# envs/monsda.yaml pins "monsda =<version>" for bioconda. Rewrite that pin to
# the version being released instead of relying on it being bumped by hand,
# so this only ever fails because bioconda does not have the release yet, not
# because the committed pin is stale.
RUN if [ "${ENVIRONMENT}" = "monsda" ]; then \
      sed -i "s/^\( *-\s*monsda\s*\)=.*/\1=${MONSDA_VERSION}/" /tmp/environment.yaml; \
    fi

RUN conda env create --yes --prefix /opt/conda/envs/monsda --file /tmp/environment.yaml \
    && conda clean --all --yes

ENV PATH="/opt/conda/envs/monsda/bin:${PATH}"
ENV MONSDA_VERSION="${MONSDA_VERSION}"

USER monsda
CMD ["/bin/bash"]
