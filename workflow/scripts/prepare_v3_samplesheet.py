#!/usr/bin/env python3
"""Build the samplesheet for the WHO v3 rerun.

Union of:
  1. Goig et al. NEJM 2025 public genomes (one row per genome ID, runs merged)
  2. Phenotyped Illumina runs (illumina.samplesheet.csv) whose run is not in (1)
  3. HDR runs (nejmc88k.csv) not already in (1) or (2)

Rows whose accession is not a run (SRX/ERS/PRJ, blank) are resolved to runs via
ENA when --resolve is given (needs internet: run on the login node). Without it,
or if ENA can't resolve them, they go to the exclusions file with a reason.

Usage:
  python workflow/scripts/prepare_v3_samplesheet.py [--resolve]
Writes:
  config/samplesheets/v3_rerun.samplesheet.csv
  config/samplesheets/v3_rerun.excluded.csv
"""
import argparse
import re
import sys
import time

import pandas as pd

GOIG = "config/samplesheets/Goig_NEJM_2025_genomes_country_accessions.tsv"
PHENO = "config/illumina.samplesheet.csv"
HDR = "config/samplesheets/nejmc88k.csv"
OUT = "config/samplesheets/v3_rerun.samplesheet.csv"
EXCL = "config/samplesheets/v3_rerun.excluded.csv"

RUN = re.compile(r"^[EDS]RR\d+$")
ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport"


def split_runs(s):
    toks = [t.strip() for t in re.split(r"[;\s]+", s if isinstance(s, str) else "") if t.strip()]
    return [t for t in toks if RUN.match(t)], [t for t in toks if not RUN.match(t)]


def ena_runs(accession):
    import requests

    for attempt in range(3):
        try:
            r = requests.get(ENA, params={"accession": accession, "result": "read_run",
                                          "fields": "run_accession", "format": "tsv"}, timeout=60)
            if r.status_code in (204, 400, 404):
                return []
            r.raise_for_status()
            return [l.split("\t")[0] for l in r.text.splitlines()[1:] if RUN.match(l.split("\t")[0])]
        except Exception as e:
            print(f"  ENA error for {accession} (attempt {attempt + 1}): {e}", file=sys.stderr)
            time.sleep(2)
    return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolve", action="store_true", help="resolve non-run accessions via ENA")
    args = ap.parse_args()

    rows, excluded = [], []

    # 1. Goig genomes
    g = pd.read_csv(GOIG, sep="\t", header=None, dtype=str,
                    names=["sample_id", "country", "bioproject", "biosample", "runs"])
    for r in g.itertuples(index=False):
        runs, other = split_runs(r.runs)
        if not runs and args.resolve:
            # prefer the biosample (sample-specific); a PRJ accession would pull a whole project
            candidates = [r.biosample] if isinstance(r.biosample, str) else []
            candidates += [t for t in other if not t.startswith("PRJ")]
            for acc in candidates:
                runs = ena_runs(acc)
                if runs:
                    break
        if not runs:
            excluded.append(dict(sample_id=r.sample_id, biosample=r.biosample, accession=r.runs,
                                 source="goig", reason="no run accession" + ("" if args.resolve else " (not resolved)")))
            continue
        rows.append(dict(sample_id=r.sample_id, biosample=r.biosample, country=r.country,
                         runs=";".join(dict.fromkeys(runs)), source="goig"))

    df = pd.DataFrame(rows)
    # identical run sets listed under two genome IDs -> profile once
    dup = df.runs.duplicated()
    for r in df[dup].itertuples(index=False):
        excluded.append(dict(sample_id=r.sample_id, biosample=r.biosample, accession=r.runs,
                             source="goig", reason="duplicate of another genome ID"))
    df = df[~dup]
    seen = {x for s in df.runs for x in s.split(";")}

    # 2. phenotyped runs not in Goig
    p = pd.read_csv(PHENO, usecols=["run", "biosample"], dtype=str)
    p = p[~p.run.isin(seen)].drop_duplicates("run")
    df = pd.concat([df, pd.DataFrame(dict(sample_id=p.run, biosample=p.biosample, country=None,
                                          runs=p.run, source="phenotyped"))])
    seen |= set(p.run)

    # 3. HDR runs not yet covered
    h = pd.read_csv(HDR, dtype=str)
    hdr_rows = []
    for r in h.itertuples(index=False):
        runs = [x for x in split_runs(r.Accession)[0]]
        if runs and not set(runs) & seen:
            hdr_rows.append(dict(sample_id=r.BioSample, biosample=r.BioSample, country=r.COUNTRY_ISOLATION,
                                 runs=";".join(runs), source="hdr"))
            seen |= set(runs)
    df = pd.concat([df, pd.DataFrame(hdr_rows)])

    assert df.sample_id.is_unique, "sample_id collision"
    df.to_csv(OUT, index=False)
    pd.DataFrame(excluded).to_csv(EXCL, index=False)
    print(df.source.value_counts().to_string())
    print(f"total samples: {len(df)}, runs: {len(seen)}; excluded: {len(excluded)} -> {EXCL}")


if __name__ == "__main__":
    main()
