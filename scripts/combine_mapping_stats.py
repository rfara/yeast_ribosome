from pathlib import Path
import csv


fieldnames = [
    "sample",
    "total_primary_reads",
    "mapped_primary_reads",
    "unmapped_primary_reads",
    "mapped_fraction",
]

Path(snakemake.output.summary).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.output.summary, "w", newline="") as out_handle:
    writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for stats_file in snakemake.input:
        with open(stats_file, newline="") as in_handle:
            reader = csv.DictReader(in_handle, delimiter="\t")
            writer.writerow(next(reader))

Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).write_text(f"merged_samples\t{len(snakemake.input)}\n")
