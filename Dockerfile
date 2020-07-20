FROM nfcore/base:1.9
LABEL authors="Marta Pozuelo" \
      description="Docker image containing all requirements for the Flomics/rnaseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-rnaseq-1.4.2/bin:$PATH


# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
