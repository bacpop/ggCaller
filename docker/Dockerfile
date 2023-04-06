FROM mambaorg/micromamba:latest

USER root

# create a project directory inside user home
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY . /src
WORKDIR /src

COPY --chown=$user:$user ../environment_linux.yml /tmp/environment_linux.yml
RUN micromamba install -y -n base -f /tmp/environment_linux.yml && \
    micromamba clean --all --yes && python -m pip install --no-deps --ignore-installed .

WORKDIR /workdir
