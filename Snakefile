from collections import defaultdict
from pathlib import Path
from itertools import product
from snakemake.utils import min_version
import pandas as pd
import glob

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

# =====================================
target_files = set()
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
