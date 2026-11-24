#!/usr/bin/env python3
"""Flag runs whose ENA FASTQ layout would trip the pipeline.
Usage: python check_ena_layout.py samplesheet.csv run_column out.tsv
"""
import sys
import io
import time
import pandas as pd
import requests

sheet, col, out = sys.argv[1:4]
runs = pd.read_csv(sheet, usecols=[col], low_memory=False)[
    col].dropna().astype(str).unique()

URL = "https://www.ebi.ac.uk/ena/portal/api/search"
FIELDS = "run_accession,library_layout,base_count,fastq_ftp,fastq_bytes"
frames = []
for i in range(0, len(runs), 500):
    batch = runs[i:i + 500]
    r = requests.post(URL, data={
        "result": "read_run",
        "includeAccessions": ",".join(batch),
        "fields": FIELDS,
        "format": "tsv",
    }, timeout=120)
    r.raise_for_status()
    frames.append(pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str))
    print(f"{min(i + 500, len(runs))}/{len(runs)}", file=sys.stderr)
    time.sleep(0.5)

df = pd.concat(frames, ignore_index=True)
df["n_files"] = df["fastq_ftp"].fillna("").apply(
    lambda s: len([x for x in s.split(";") if x]))
has_pair = df["fastq_ftp"].fillna("").str.contains(
    r"_1\.fastq\.gz") & df["fastq_ftp"].fillna("").str.contains(r"_2\.fastq\.gz")

df["issue"] = ""
df.loc[(df.n_files == 3), "issue"] = "paired+orphan (3 files)"
df.loc[(df.library_layout == "PAIRED") & ~has_pair,
       "issue"] = "paired but no _1/_2 on ENA"
df.loc[df.n_files == 0, "issue"] = "no fastq on ENA"
missing = set(runs) - set(df.run_accession)

flagged = df[df.issue != ""]
flagged.to_csv(out, sep="\t", index=False)
print(f"{len(runs)} runs queried, {len(missing)} not returned by ENA")
print(flagged.issue.value_counts().to_string())
