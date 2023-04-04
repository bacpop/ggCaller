FROM mambaorg/micromamba:latest

ARG USER_ID
ARG GROUP_ID

USER root

RUN addgroup --gid $GROUP_ID user
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user
USER user

# create a project directory inside user home
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY . /src
WORKDIR /src

COPY --chown=$user:$user ../environment_linux.yml /tmp/environment_linux.yml
RUN micromamba install -y -n base -f /tmp/environment_linux.yml && \
    micromamba clean --all --yes && python -m pip install --no-deps --ignore-installed .

WORKDIR /workdir
