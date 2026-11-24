#!/usr/bin/env python3
"""Drop non-WGS runs (e.g. AMPLICON, Targeted-Capture) from the v3 samplesheet.

Amplicon/targeted runs only cover a few genes, so TB-Profiler reports every other
drug as susceptible and the genome median depth is ~0. Profiling them wastes jobs
and, if filtered on target depth instead of genome depth, gives false negatives.

Queries ENA for library_strategy/library_source of every run (needs internet:
run on the login node). Within a sample, non-WGS runs are dropped; a sample with
no WGS runs left is moved to the exclusions file. Runs ENA doesn't return are kept.

Usage: python workflow/scripts/filter_wgs_runs.py
Rewrites config/samplesheets/v3_rerun.samplesheet.csv (backup: .pre_wgs.csv),
appends to v3_rerun.excluded.csv and writes v3_rerun.run_metadata.tsv.
"""
import io
import shutil
import sys
import time

import pandas as pd
import requests

SHEET = "config/samplesheets/v3_rerun.samplesheet.csv"
EXCL = "config/samplesheets/v3_rerun.excluded.csv"
META = "config/samplesheets/v3_rerun.run_metadata.tsv"
URL = "https://www.ebi.ac.uk/ena/portal/api/search"
FIELDS = "run_accession,library_strategy,library_source,library_layout,instrument_platform,base_count"
BATCH = 500

df = pd.read_csv(SHEET, dtype=str)
runs = sorted({r for s in df.runs for r in s.split(";")})

frames = []
for i in range(0, len(runs), BATCH):
    batch = runs[i:i + BATCH]
    for attempt in range(5):
        try:
            r = requests.post(URL, data={"result": "read_run", "format": "tsv", "fields": FIELDS,
                                         "includeAccessions": ",".join(batch)}, timeout=180)
            r.raise_for_status()
            break
        except Exception as e:
            print(f"batch {i // BATCH}: {e}; retrying", file=sys.stderr)
            time.sleep(10 * (attempt + 1))
    else:
        sys.exit(f"ENA query failed for batch starting {batch[0]}")
    if r.text.strip():
        frames.append(pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str))
    print(f"{min(i + BATCH, len(runs))}/{len(runs)} runs queried", file=sys.stderr)
    time.sleep(0.3)

meta = pd.concat(frames, ignore_index=True).drop_duplicates("run_accession")
meta.to_csv(META, sep="\t", index=False)
strategy = dict(zip(meta.run_accession, meta.library_strategy.fillna("")))
not_found = [r for r in runs if r not in strategy]

keep_rows, excluded = [], []
for row in df.itertuples(index=False):
    rs = row.runs.split(";")
    wgs = [r for r in rs if strategy.get(r, "WGS") == "WGS"]   # unknown -> keep
    if wgs:
        keep_rows.append(row._replace(runs=";".join(wgs)))
    else:
        kinds = ",".join(sorted({strategy[r] for r in rs}))
        excluded.append(dict(sample_id=row.sample_id, biosample=row.biosample, accession=row.runs,
                             source=row.source, reason=f"not WGS ({kinds})"))

shutil.copy(SHEET, SHEET.replace(".csv", ".pre_wgs.csv"))
pd.DataFrame(keep_rows, columns=df.columns).to_csv(SHEET, index=False)
old = pd.read_csv(EXCL, dtype=str) if pd.io.common.file_exists(EXCL) else pd.DataFrame()
(pd.concat([old, pd.DataFrame(excluded)], ignore_index=True)
   .drop_duplicates(["sample_id", "reason"]).to_csv(EXCL, index=False))

print(meta.library_strategy.value_counts(dropna=False).to_string())
print(f"\nruns not returned by ENA (kept): {len(not_found)}")
print(f"samples kept: {len(keep_rows)} / {len(df)}; excluded as not WGS: {len(excluded)}")
