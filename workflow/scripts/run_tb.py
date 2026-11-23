import pandas as pd
from pathlib import Path 
import subprocess 
import tempfile

run_p = Path('/data/scratch/projects/punim1703/false_negative/old_tbprof/resources/data/samples.csv')
samples = pd.read_csv(run_p)
runs = samples['run'].unique()

for run in runs:
    out_dir = Path(f'/data/scratch/projects/punim1703/false_negative/tbprof-rerun/results/runs/{run}')
    res_jsn = next((out_dir/"results").glob("*.json"), None)
    res_txt = next((out_dir/"results").glob("*.txt"), None)
    vcf = next((out_dir/"vcf").glob("*.targets.vcf.gz"), None)

    if (
        out_dir.exists()
        and res_jsn is not None
        and res_txt is not None 
        and vcf is not None 
    ):
        continue
    
    out_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:

        try:
            subprocess.run(
            ["fastq-dl", "-a", run, "--outdir", str(tmpdir)],
            check=True
            )
        except subprocess.CalledProcessError as e:
            print(e)
            continue

        
        reads = list(Path(tmpdir).glob("*.gz"))
        reads = [str(f) for f in reads] 

        if not reads:
            print('reads not downloaded')
            continue

        if len(reads) < 2:
            interleaved = reads[0]
            r1 = str(Path(tmpdir) / f"{run}_R1.fastq.gz")
            r2 = str(Path(tmpdir) / f"{run}_R2.fastq.gz")

            try:
                subprocess.run(
                ["seqfu", "deinterleave", interleaved, "-o", r1, r2],
                check=True
            )
            except subprocess.CalledProcessError as e:
                print(e)
                continue

            reads = [r1, r2]

        try:
            subprocess.run([
                        'tb-profiler', 'profile',
                        '-1', f"{reads[0]}",
                        '-2', f"{reads[1]}",
                        '--txt',
                        '--no_trim',
                        '-t', '12',
                        '-d', str(out_dir)
                    ], check=True)
        except subprocess.CalledProcessError as e:
            print(e)
            continue
        
        try:
            subprocess.run([
                        'rm',
                        '-rf',
                        f'{out_dir}/bam'
                    ], check=True)
        except subprocess.CalledProcessError as e:
            print(e)
            continue