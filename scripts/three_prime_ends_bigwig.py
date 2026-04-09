import numpy as np
import pysam
import pyBigWig


def main(input_bam, chrom_sizes_fai, output_bw, log_file):
    chroms = {}
    with open(chrom_sizes_fai) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            chroms[parts[0]] = int(parts[1])

    counts = {chrom: np.zeros(length, dtype=np.float64) for chrom, length in chroms.items()}
    total_reads = 0

    with pysam.AlignmentFile(input_bam, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            end_pos = read.reference_end - 1
            chrom = read.reference_name
            if chrom in counts and 0 <= end_pos < len(counts[chrom]):
                counts[chrom][end_pos] += 1
                total_reads += 1

    bw = pyBigWig.open(output_bw, "w")
    bw.addHeader(list(chroms.items()))

    for chrom, length in chroms.items():
        arr = counts[chrom]
        nonzero = np.nonzero(arr)[0]
        if len(nonzero) == 0:
            continue
        starts = nonzero.astype(np.int64).tolist()
        ends = (nonzero + 1).astype(np.int64).tolist()
        values = arr[nonzero].tolist()
        bw.addEntries([chrom] * len(starts), starts, ends=ends, values=values)

    bw.close()

    with open(log_file, "w") as lf:
        lf.write(f"total_reads_counted\t{total_reads}\n")
        for chrom in chroms:
            lf.write(f"{chrom}_total_signal\t{counts[chrom].sum()}\n")


if __name__ == "__main__":
    main(
        input_bam=snakemake.input.bam,
        chrom_sizes_fai=snakemake.input.fai,
        output_bw=snakemake.output.bw,
        log_file=snakemake.log[0],
    )
