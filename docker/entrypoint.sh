#!/bin/bash --login
set -e

micromamba activate base
exec "$@"