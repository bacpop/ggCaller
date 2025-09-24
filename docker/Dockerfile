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
RUN apt-get update && apt-get install -y build-essential wget tar && \
    wget http://ftp.debian.org/debian/pool/main/e/execstack/execstack_1.0.0.orig.tar.gz && \
    tar xf execstack_1.0.0.orig.tar.gz && \
    cd execstack-1.0.0 && make && cp execstack /usr/local/bin/ && cd .. && rm -rf execstack-1.0.0* && \
    rm -rf /var/lib/apt/lists/*

# Clear executable stack flag for libtorch_cpu.so (PyTorch 1.10)
RUN execstack -c $(python -c "import torch, os; print(os.path.join(os.path.dirname(torch.__file__), 'lib', 'libtorch_cpu.so'))")

# download ggCaller DB during build
RUN ggcaller --balrog-db ggc_db

ENV PATH=$PATH:/opt/conda/bin

ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]
