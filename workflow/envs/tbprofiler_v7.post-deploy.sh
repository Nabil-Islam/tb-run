#!/usr/bin/env bash
# Run by Snakemake (--use-conda) once, after it creates the tbprofiler_v7 env.
# Installs TB-Profiler 7.0.0 and its default WHO v3 database, mirroring the
# bioconda recipe's build.sh.
#
# Every tool is called by its full path inside the env, so this works even if
# Snakemake activates the env with a different conda (e.g. a cluster module).
set -euo pipefail

# Snakemake copies this script to "<env prefix>.post-deploy.sh"
ENV="$(cd "$(dirname "$0")" && pwd)/$(basename "$0" .post-deploy.sh)"
if [ ! -x "$ENV/bin/python" ]; then
    echo "post-deploy: cannot find the env python at $ENV/bin/python" >&2
    exit 1
fi
export PATH="$ENV/bin:$PATH"
export PYTHONNOUSERSITE=1   # never install into or import from ~/.local

VERSION="7.0.0"
src=$(mktemp -d)
trap 'rm -rf "$src"' EXIT

git -c advice.detachedHead=false clone --quiet --depth 1 --branch "v${VERSION}" \
    https://github.com/jodyphelan/TBProfiler "$src/TBProfiler"
cd "$src/TBProfiler"
"$ENV/bin/python" -m pip install --no-deps --ignore-installed --no-build-isolation --no-user .

cd db
"$ENV/bin/samtools" dict who_v3/genome.fasta -o who_v3/genome.dict
"$ENV/bin/samtools" faidx who_v3/genome.fasta
"$ENV/bin/bwa" index who_v3/genome.fasta
"$ENV/bin/python" "$ENV/bin/tb-profiler" load_library who_v3 --force

"$ENV/bin/python" "$ENV/bin/tb-profiler" version
"$ENV/bin/python" "$ENV/bin/tb-profiler" list_db
