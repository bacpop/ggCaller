#!/bin/bash --login
set -e

micromamba activate $HOME/app/env
exec "$@"