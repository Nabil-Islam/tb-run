#!/usr/bin/env bash
# Run by Snakemake (--use-conda) inside the new tbprofiler_v7 env, once, after it
# is created. Installs TB-Profiler 7.0.0 and its default WHO v3 database, mirroring
# the bioconda recipe's build.sh.
set -euo pipefail

VERSION="7.0.0"
src=$(mktemp -d)
trap 'rm -rf "$src"' EXIT

git clone --quiet --depth 1 --branch "v${VERSION}" https://github.com/jodyphelan/TBProfiler "$src/TBProfiler"
cd "$src/TBProfiler"
python -m pip install --no-deps --ignore-installed --no-build-isolation .

cd db
samtools dict who_v3/genome.fasta -o who_v3/genome.dict
samtools faidx who_v3/genome.fasta
bwa index who_v3/genome.fasta
tb-profiler load_library who_v3 --force

tb-profiler version
tb-profiler list_db
