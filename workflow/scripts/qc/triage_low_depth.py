#!/usr/bin/env python3
"""Usage: python triage_low_depth.py low_depth_all.txt triage.tsv"""
import sys
import io
import requests
import pandas as pd

src, out = sys.argv[1:3]
rows = [l.split() for l in open(src).read().splitlines()[1:] if l.strip()]
df = pd.DataFrame(rows, columns=["run", "pct_mapped", "depth"])
df[["pct_mapped", "depth"]] = df[["pct_mapped", "depth"]].apply(
    pd.to_numeric, errors="coerce")

ena = []
for i in range(0, len(df), 500):
    r = requests.post("https://www.ebi.ac.uk/ena/portal/api/search", data={
        "result": "read_run", "format": "tsv",
        "includeAccessions": ",".join(df.run[i:i+500]),
        "fields": "run_accession,study_accession,library_layout,base_count,fastq_ftp",
    }, timeout=120)
    r.raise_for_status()
    ena.append(pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str))
ena = pd.concat(ena).rename(columns={"run_accession": "run"})
df = df.merge(ena, on="run", how="left")
df["expected_depth"] = df.base_count.astype(float) / 4_411_532


def cat(r):
    if r.pct_mapped == 0:
        return "zero_mapped_RERUN" if r.expected_depth >= 10 else "zero_mapped_tiny_submission"
    if r.pct_mapped < 80:
        return "low_mapping_contaminated"
    return "low_sequencing_depth"


df["category"] = df.apply(cat, axis=1)
df.to_csv(out, sep="\t", index=False)
print(df.category.value_counts().to_string())
print("\nzero_mapped_RERUN by study:")
print(df[df.category == "zero_mapped_RERUN"].study_accession.value_counts().head(
    15).to_string())
