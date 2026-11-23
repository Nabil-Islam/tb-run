"""Combine TB-Profiler 7 (WHO v3) JSON reports into one long CSV.

Same core columns as combine_tbprofiler_reports.py (run, biosample, bioproject,
technology, tool, drug, prediction, lineage, mutations), so downstream code keeps
working, plus QC and database columns.

Calls come from `drug_resistance_interpretation` (v7), which includes lineage-based
intrinsic resistance and marks drugs whose only variants carry warnings as "W".
Older JSONs without that field fall back to `dr_variants`.
"""
import sys

sys.stderr = open(snakemake.log[0], "w")

import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

QC_FIELDS = [
    "num_reads_mapped",
    "percent_reads_mapped",
    "genome_median_depth",
    "target_median_depth",
]
HEADER = [
    "run", "biosample", "bioproject", "technology", "tool", "drug",
    "prediction", "lineage", "mutations", "resistance_level",
    *QC_FIELDS, "tbprofiler_version", "db_name", "db_commit", "sample_id",
]


def calls_from_interpretation(data):
    """drug -> (prediction, mutations, level) for drugs with any call."""
    out = {}
    for d in data["drug_resistance_interpretation"]:
        interp = d.get("interpretation", "")
        if not interp:
            continue
        pred = "W" if interp == "W" else "R"
        level = ""
        if interp.startswith("R ("):
            level = interp[3:-1]  # "low" / "high"
        # "rpoB:p.Ser450Leu" -> "rpoB_p.Ser450Leu", matching the v6 summary format
        muts = ";".join(sorted({m.replace(":", "_", 1) for m in d.get("resistance_mechanisms", [])}))
        out[d["drug"]] = (pred, muts, level)
    return out


def calls_from_dr_variants(data):
    drug_variants = defaultdict(set)
    for v in data.get("dr_variants", []):
        for info in v.get("drugs", []):
            if info.get("type", "drug_resistance") == "drug_resistance":
                drug_variants[info["drug"]].add(f"{v['gene_name']}_{v['change']}")
    return {d: ("R", ";".join(sorted(m)), "") for d, m in drug_variants.items()}


runs_by_sample = snakemake.params.runs
biosample_by_sample = getattr(snakemake.params, "biosamples", None) or {}
tech = snakemake.params.tech
n = 0

with open(snakemake.output.report, "w", newline="") as fout:
    w = csv.writer(fout)
    w.writerow(HEADER)
    for p in map(Path, snakemake.input.reports):
        sample = p.parent.name
        opener = gzip.open if p.suffix == ".gz" else open
        with opener(p, "rt") as fh:
            data = json.load(fh)

        if "drug_resistance_interpretation" in data:
            calls = calls_from_interpretation(data)
        else:
            calls = calls_from_dr_variants(data)

        qc = data.get("qc", {}) or {}
        pipeline = data.get("pipeline", {}) or {}
        db = pipeline.get("db_version", data.get("db_version", {})) or {}
        common = dict(
            run=runs_by_sample.get(sample, ""),
            biosample=biosample_by_sample.get(sample) or sample,
            sample_id=sample,
            bioproject="",
            technology=tech,
            tool="tbprofiler",
            lineage=data.get("sub_lineage", ""),
            tbprofiler_version=pipeline.get("software_version", data.get("tbprofiler_version", "")),
            db_name=db.get("name", ""),
            db_commit=db.get("commit", ""),
            **{k: qc.get(k, "") for k in QC_FIELDS},
        )

        if not calls:
            rows = [dict(common, drug="all", prediction="S", mutations="", resistance_level="")]
        else:
            rows = [
                dict(common, drug=drug, prediction=pred, mutations=muts, resistance_level=lvl)
                for drug, (pred, muts, lvl) in sorted(calls.items())
            ]
        for r in rows:
            w.writerow([r[h] for h in HEADER])
        n += 1

print(f"combined {n} reports", file=sys.stderr)
