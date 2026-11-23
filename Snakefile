from collections import defaultdict
from pathlib import Path
from itertools import product
from snakemake.utils import min_version
import pandas as pd
import glob
import sys

# =====================================
# Workflow config items
configfile: "config/config.yaml"

WORKFLOW = Path("workflow").resolve()
RULES = WORKFLOW / "rules"
ENVS = WORKFLOW / "envs"
SCRIPTS = WORKFLOW / "scripts"
RESULTS = Path("results").resolve()
PLOTS = RESULTS / "plots"
TABLES = RESULTS / "tables"
CAPTIONS = WORKFLOW / "report"
RESOURCES = Path("resources").resolve()
LOGS = Path("logs/rules").resolve()
BENCH = RESULTS / "benchmark"
CONTAINERS = config["containers"]
GB = 1_024
MIN_ILLUMINA_COV = config["min_cov"]
MIN_NANOPORE_COV = config["min_cov"]

TECHS = ["nanopore", "illumina"]
TOOLS = ["tbprofiler"]

inclusion_expr = f"illumina_covg >= {MIN_ILLUMINA_COV} and nanopore_covg >= {MIN_NANOPORE_COV} and lineage != 'mixed'"
h2h_df = pd.read_csv(config["h2h_samplesheet"]).query(inclusion_expr)

illumina_df = pd.read_csv(
    config["illumina_samplesheet"], low_memory=False, index_col="run"
)
illumina_df["run"] = illumina_df.index
# illumina_query = "run in @illumina_subsample_runs or delamanid.notna()"
# illumina_df.query(illumina_query, inplace=True)
illumina_depth_df = illumina_df.query("(levofloxacin.notna() and streptomycin.notna() and kanamycin.isna() and ofloxacin.notna() and pyrazinamide.notna()) or (delamanid.notna() and pyrazinamide.notna() and levofloxacin.isna())")
ont_df = pd.read_csv(config["nanopore_samplesheet"], index_col="run", low_memory=False)
# genomes from 88k genomes yet to be profiled 
ext_df = pd.read_csv(config["extended_samplesheet"], low_memory=False, index_col="Accession")
ext_df["Accession"] = ext_df.index
EXT_PROJ = "extended"   # no BioProject column, so use a fixed folder name
EXT_TECH = "illumina"   # tb-profiler --platform; must be illumina or nanopore

# stop {tech} from matching anything else
wildcard_constraints:
    tech="illumina|nanopore",
    sample="[^/]+",

# =====================================
# WHO v3 rerun (TB-Profiler 7.0.0): one job per sample, runs merged
V3 = config["v3_rerun"]
V3_DIR = RESULTS / f"amr_predictions/tbprofiler_{V3['tag']}" / V3["tech"]
V3_RUNS = {}
if Path(V3["samplesheet"]).exists():
    _v3 = pd.read_csv(V3["samplesheet"], low_memory=False, dtype=str)
    _v3 = _v3.dropna(subset=[V3["sample_col"], V3["runs_col"]])
    _dups = _v3[V3["sample_col"]][_v3[V3["sample_col"]].duplicated()].unique()
    if len(_dups):
        raise ValueError(f"{len(_dups)} duplicated sample IDs in {V3['samplesheet']}, e.g. {list(_dups[:5])}")
    for _s, _r in zip(_v3[V3["sample_col"]], _v3[V3["runs_col"]]):
        V3_RUNS[_s.strip()] = [x.strip() for x in _r.split(";") if x.strip()]
else:
    print(f"WARNING: {V3['samplesheet']} not found; WHO v3 rerun targets skipped", file=sys.stderr)

# =====================================
target_files = set()
if V3_RUNS:
    target_files.add(V3_DIR / "summary.csv")
if config.get("run_v6", False):   # original TB-Profiler 6.6.6 runs
    for tech in TECHS:
        target_files.add(RESULTS / f"amr_predictions/tbprofiler/{tech}/summary.csv")

# =====================================
rule all:
    input:
        target_files, 

# =====================================

rule download_verify_predict:
    output:
        report=RESULTS / "amr_predictions/tbprofiler/{tech}/{proj}/{sample}/{run}/results/{run}.results.json",
        run_info=RESULTS / "validate/{tech}/{proj}/{sample}/{run}/run_info.tsv",
    log:
        LOGS / "download_verify_predict/{tech}/{proj}/{sample}/{run}/{run}.log",
    conda:
        ENVS/'tbprofiler.yaml'
    threads: 4 
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(8 * GB),
        runtime="4h"
    params:
        db="ena",
        fq_opts="--verbose",
        tb_opts="--txt --no_trim -p {run}",
        outdir=lambda wildcards, output: Path(output.report).parent.parent,
        script_dir=SCRIPTS,
        res_dir=RESULTS
    group:
        "qc"
    shell:
        """
        tmpdir=$(mktemp -d)
        trap 'rm -rf "$tmpdir"' EXIT

        cd "$tmpdir"

        shopt -s nullglob

        clean_bad_reads() {{
            for f in *.fastq.gz; do
                gzip -t "$f" 2>/dev/null || {{ echo "removing corrupt $f" >> {log}; rm -f "$f"; }}
            done
        }}

        fastq-dl {params.fq_opts} -a {wildcards.run} --provider {params.db} 2>> {log} || true
        clean_bad_reads
        reads=(*.fastq.gz)

        if [ ${{#reads[@]}} -lt 1 ]; then
            fastq-dl {params.fq_opts} -a {wildcards.run} --provider sra 2>> {log} || true
            clean_bad_reads
            reads=(*.fastq.gz)
        fi

        if [ ${{#reads[@]}} -lt 1 ]; then
            echo "No reads downloaded for {wildcards.run}" >> {log}
            exit 1
        fi

        python3 {params.script_dir}/validate_run_info.py fastq-run-info.tsv {wildcards.run} {wildcards.tech} {params.res_dir}

        exec 2>> {log} 

        run_acc={wildcards.run}
        run_info=run_info.tsv 

        r1={wildcards.run}_1.fastq.gz
        r2={wildcards.run}_2.fastq.gz
        se={wildcards.run}.fastq.gz

        if [ -s "$r1" ] && [ -s "$r2" ]; then
            input_arg=("-1" "$r1" "-2" "$r2")      # paired; ignore any orphan file
        elif [ -s "$se" ]; then
            input_arg=("-1" "$se")                 # genuine single-end (or nanopore)
        else
            echo "No usable reads for {wildcards.run}: $(ls)" >> {log}
            exit 1
        fi

        tb-profiler profile "${{input_arg[@]}}" {params.tb_opts} -t {threads} --platform {wildcards.tech}
        
        mkdir -p $(dirname {output.report})
        mkdir -p $(dirname {output.run_info})

        mv $tmpdir/results/{wildcards.run}.results.json {output.report}
        mv run_info.tsv {output.run_info}

        """


def infer_tbprofiler_reports(wildcards):
    if wildcards.tech == "illumina":
        df = illumina_df
    else:
        df = ont_df

    files = []
    for run, row in df.iterrows():
        proj = row["bioproject"]
        sample = row["biosample"]
        p = (
            RESULTS
            / f"amr_predictions/tbprofiler/{wildcards.tech}/{proj}/{sample}/{run}/results/{run}.results.json"
        )
        files.append(p)

    return files

def infer_ext_tbprofiler_reports(wildcards):
    return [
        RESULTS
        / f"amr_predictions/tbprofiler/{EXT_TECH}/{EXT_PROJ}/{row['BioSample']}/{run}/results/{run}.results.json"
        for run, row in ext_df.iterrows()
    ]

rule combine_ext_tbprofiler_reports:
    input:
        reports=infer_ext_tbprofiler_reports,
    output:
        report=RESULTS / f"amr_predictions/tbprofiler/{EXT_TECH}/extended_summary.csv",
    log:
        LOGS / "combine_tbprofiler_reports/extended.log",
    resources:
        mem_mb=GB,
    script:
        str(SCRIPTS / "combine_tbprofiler_reports.py")

rule combine_tbprofiler_reports:
    input:
        reports=infer_tbprofiler_reports,
    output:
        report=RESULTS / "amr_predictions/tbprofiler/{tech}/summary.csv",
    log:
        LOGS / "combine_tbprofiler_reports/{tech}.log",
    resources:
        mem_mb=GB,
    script:
        str(SCRIPTS / "combine_tbprofiler_reports.py")


# =====================================
# WHO v3 rerun

rule profile_sample_v3:
    output:
        report=V3_DIR / "{sample}/{sample}.results.json",
        vcf=V3_DIR / "{sample}/{sample}.targets.vcf.gz",
        run_info=V3_DIR / "{sample}/run_info.tsv",
    log:
        LOGS / f"profile_sample_v3/{V3['tag']}/{{sample}}.log",
    conda:
        ENVS / "tbprofiler_v7.yaml"
    threads: 4
    retries: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(8 * GB),
        runtime="4h",
    params:
        runs=lambda wildcards: " ".join(V3_RUNS[wildcards.sample]),
        db=V3["db"],
        tech=V3["tech"],
    group:
        "qc"
    shell:
        """
        exec 2>> {log}
        echo "== $(date) sample={wildcards.sample} runs={params.runs}" >&2
        tmpdir=$(mktemp -d -p "${{TMPDIR:-/tmp}}")
        trap 'rm -rf "$tmpdir"' EXIT
        cd "$tmpdir"
        shopt -s nullglob

        r1s=(); r2s=(); ses=()
        for run in {params.runs}; do
            d="dl/$run"; mkdir -p "$d"
            for provider in ena sra; do
                fastq-dl -a "$run" --provider $provider --only-provider -o "$d" --verbose >&2 || true
                for f in "$d"/*.fastq.gz; do
                    gzip -t "$f" 2>/dev/null || {{ echo "removing corrupt $f" >&2; rm -f "$f"; }}
                done
                if {{ [ -s "$d/${{run}}_1.fastq.gz" ] && [ -s "$d/${{run}}_2.fastq.gz" ]; }} || [ -s "$d/${{run}}.fastq.gz" ]; then
                    break
                fi
                echo "no usable reads for $run from $provider" >&2
                rm -f "$d"/*.fastq.gz
            done

            # run metadata (keep one header) + platform check
            if [ -s "$d/fastq-run-info.tsv" ]; then
                [ -s run_info.tsv ] || head -n 1 "$d/fastq-run-info.tsv" > run_info.tsv
                tail -n +2 "$d/fastq-run-info.tsv" >> run_info.tsv
                platform=$(python3 -c "import csv,sys; r=list(csv.DictReader(open(sys.argv[1]),delimiter='\\t')); print(r[0].get('instrument_platform',''))" "$d/fastq-run-info.tsv")
                case "${{platform^^}}" in
                    *ILLUMINA*) expected=illumina ;;
                    *OXFORD*|*NANOPORE*) expected=nanopore ;;
                    *) expected=unknown ;;
                esac
                if [ "$expected" != "{params.tech}" ]; then
                    echo "platform '$platform' for $run does not match {params.tech}" >&2
                    exit 1
                fi
            fi

            if [ -s "$d/${{run}}_1.fastq.gz" ] && [ -s "$d/${{run}}_2.fastq.gz" ]; then
                r1s+=("$d/${{run}}_1.fastq.gz"); r2s+=("$d/${{run}}_2.fastq.gz")
            elif [ -s "$d/${{run}}.fastq.gz" ]; then
                ses+=("$d/${{run}}.fastq.gz")
            else
                echo "No usable reads for $run: $(ls "$d")" >&2
                exit 1
            fi
        done

        if [ ${{#r1s[@]}} -gt 0 ] && [ ${{#ses[@]}} -gt 0 ]; then
            echo "mixed paired/single-end runs for {wildcards.sample}; not merging" >&2
            exit 1
        fi

        # merge runs (concatenated gzip is valid gzip)
        if [ ${{#r1s[@]}} -gt 0 ]; then
            cat "${{r1s[@]}}" > R1.fastq.gz; cat "${{r2s[@]}}" > R2.fastq.gz
            rm -f "${{r1s[@]}}" "${{r2s[@]}}"
            input_arg=("-1" "R1.fastq.gz" "-2" "R2.fastq.gz")
        else
            cat "${{ses[@]}}" > SE.fastq.gz; rm -f "${{ses[@]}}"
            input_arg=("-1" "SE.fastq.gz")
        fi

        tb-profiler profile "${{input_arg[@]}}" --db {params.db} --platform {params.tech} \
            -p {wildcards.sample} -t {threads} --no_trim --temp "$tmpdir" --dir "$tmpdir/out" >&2

        json="$tmpdir/out/results/{wildcards.sample}.results.json"
        mapped=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['qc'].get('num_reads_mapped', 1))" "$json")
        if [ "$mapped" = "0" ]; then
            echo "0 reads mapped for {wildcards.sample}; treating as a failed job" >&2
            exit 1
        fi

        mkdir -p "$(dirname {output.report})"
        mv "$json" {output.report}
        mv "$tmpdir/out/vcf/{wildcards.sample}.targets.vcf.gz" {output.vcf}
        if [ -s run_info.tsv ]; then mv run_info.tsv {output.run_info}; else touch {output.run_info}; fi
        """


rule combine_v3_reports:
    input:
        reports=lambda wildcards: [V3_DIR / f"{s}/{s}.results.json" for s in V3_RUNS],
    output:
        report=V3_DIR / "summary.csv",
    log:
        LOGS / f"combine_v3_reports/{V3['tag']}.log",
    params:
        runs={s: ";".join(r) for s, r in V3_RUNS.items()},
        tech=V3["tech"],
    resources:
        mem_mb=4 * GB,
        runtime="4h",
    script:
        str(SCRIPTS / "combine_v3_reports.py")
