from pathlib import Path
import subprocess


def count_primary_reads(path):
    command = ["samtools", "view", "-c", "-F", "2304", str(path)]
    return int(subprocess.check_output(command, text=True).strip())


mapped = count_primary_reads(snakemake.input.aligned)
unmapped = count_primary_reads(snakemake.input.unaligned)
total = mapped + unmapped
mapped_fraction = mapped / total if total else 0.0

Path(snakemake.output.stats).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.output.stats).write_text(
    "sample\ttotal_primary_reads\tmapped_primary_reads\tunmapped_primary_reads\tmapped_fraction\n"
    f"{snakemake.wildcards.sample}\t{total}\t{mapped}\t{unmapped}\t{mapped_fraction:.6f}\n"
)
Path(snakemake.log[0]).write_text(
    f"sample\t{snakemake.wildcards.sample}\n"
    f"total_primary_reads\t{total}\n"
    f"mapped_primary_reads\t{mapped}\n"
    f"unmapped_primary_reads\t{unmapped}\n"
    f"mapped_fraction\t{mapped_fraction:.6f}\n"
)
