#!/usr/bin/env python3
"""Usage: python check_orphan_results.py illumina_flagged.tsv results/amr_predictions/tbprofiler out.tsv"""
import sys
import json
import glob
import os
import pandas as pd

flagged_tsv, res_root, out = sys.argv[1:4]
GENOME = 4_411_532

fl = pd.read_csv(flagged_tsv, sep="\t", dtype=str)
fl = fl[fl.issue == "paired+orphan (3 files)"].copy()
fl["expected_depth"] = fl.base_count.astype(float) / GENOME

# index all result files once
paths = {os.path.basename(p).replace(".results.json", ""): p
         for p in glob.glob(f"{res_root}/*/*/*/*/results/*.results.json")}


def median_depth(path):
    qc = json.load(open(path)).get("qc", {})
    for k in ("target_median_depth", "median_depth", "region_median_depth"):
        if k in qc:
            return float(qc[k])
    return None


fl["result_exists"] = fl.run_accession.isin(paths)
fl["observed_depth"] = [median_depth(
    paths[r]) if r in paths else None for r in fl.run_accession]
fl["ratio"] = fl.observed_depth / fl.expected_depth

fl.to_csv(out, sep="\t", index=False)
print("results present:", fl.result_exists.sum(), "/", len(fl))
print(fl.loc[fl.result_exists, "ratio"].describe())
print("ratio < 0.3:", (fl.ratio < 0.3).sum())
