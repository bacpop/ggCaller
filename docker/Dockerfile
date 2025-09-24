FROM mambaorg/micromamba:latest

USER root

# create a project directory inside user home
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# project dir
ENV PROJECT_DIR=/app
RUN mkdir -p $PROJECT_DIR
COPY . $PROJECT_DIR
WORKDIR $PROJECT_DIR

# build conda env
ENV ENV_PREFIX=$PROJECT_DIR/env
COPY --chown=$user:$user docker/environment_docker.yml /tmp/environment_docker.yml

COPY --chown=$user:$user docker/entrypoint.sh /usr/local/bin/
RUN chmod u+x /usr/local/bin/entrypoint.sh

# install environment + ggCaller
RUN micromamba install -y -n base -f /tmp/environment_docker.yml && \
    micromamba clean --all --yes && \
    python -m pip install --no-deps --ignore-installed . && \
    mkdir -p ggc_db

# Build execstack from source (needed for PyTorch 1.10)
# Install patchelf to fix libtorch_cpu.so
RUN apt-get update && apt-get install -y patchelf && rm -rf /var/lib/apt/lists/*
RUN find $(python -c "import site; print(site.getsitepackages()[0])") -name libtorch_cpu.so \
    -exec patchelf --clear-execstack {} \;

# download ggCaller DB during build
RUN ggcaller --balrog-db ggc_db

ENV PATH=$PATH:/opt/conda/bin

ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]
