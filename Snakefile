from pathlib import Path


WORKFLOW_DIR = Path(workflow.basedir).resolve()
REPO_ROOT = WORKFLOW_DIR
RESULTS_DIR = WORKFLOW_DIR / "results"
LOG_DIR = WORKFLOW_DIR / "logs"
BENCHMARK_DIR = WORKFLOW_DIR / "benchmarks"
ENV_YAML = WORKFLOW_DIR / "envs" / "yeast_ribosome.yaml"
SCRIPTS_DIR = WORKFLOW_DIR / "scripts"

configfile: str(WORKFLOW_DIR / "config.yaml")

shell.executable("bash")


def repo_path(path_string):
    path = Path(path_string)
    if path.is_absolute():
        return str(path)
    return str((REPO_ROOT / path).resolve())


def executable_path(path_string):
    path = Path(path_string)
    if path.is_absolute():
        return str(path)
    if "/" in path_string:
        return str((REPO_ROOT / path).resolve())
    return path_string


def results_path(*parts):
    return str((RESULTS_DIR.joinpath(*parts)).resolve())


def log_path(*parts):
    return str((LOG_DIR.joinpath(*parts)).resolve())


def benchmark_path(*parts):
    return str((BENCHMARK_DIR.joinpath(*parts)).resolve())


def sample_fastq(wildcards):
    return repo_path(config["samples"][wildcards.sample])


def slurm_cluster_args(qos, mem, job_name, slurm_log_file, threads):
    partition, tier = qos.split("_", 1)
    time_limit = config["slurm"]["time_limits"][tier]
    parts = [
        f"--partition={partition}",
        f"--qos={qos}",
        f"--time={time_limit}",
        f"--mem={mem}",
        f"--cpus-per-task={threads}",
        f"-J {job_name}",
        f"-o {slurm_log_file}",
    ]
    extra = config["slurm"].get("sbatch_extra", "").strip()
    if extra:
        parts.append(extra)
    return " ".join(parts)


def rule_cluster_args(rule_key, job_name, slurm_log_file, threads):
    rule_cfg = config["slurm"]["jobs"][rule_key]
    return slurm_cluster_args(rule_cfg["qos"], rule_cfg["mem"], job_name, slurm_log_file, threads)


def minimap2_extra():
    parts = ["--secondary=no"]
    extra = config["minimap2"].get("extra", "").strip()
    if extra:
        parts.append(extra)
    return " ".join(parts)


SAMPLES = list(config["samples"].keys())
REFERENCE_FASTA = repo_path(config["reference"]["fasta"])
MINIMAP2_BIN = executable_path(config["minimap2"]["executable"])
MINIMAP2_PRESET = config["minimap2"]["preset"]
MINIMAP2_EXTRA = minimap2_extra()
REFERENCE_MMI = results_path("reference", "rrna_35s.mmi")
ALIGNED_BAM = results_path("aligned", "{sample}.rrna_35s.aligned.bam")
ALIGNED_BAI = results_path("aligned", "{sample}.rrna_35s.aligned.bam.bai")
UNALIGNED_BAM = results_path("unaligned", "{sample}.rrna_35s.unaligned.bam")
MAPPING_STATS = results_path("stats", "{sample}.mapping.tsv")
COMBINED_STATS = results_path("stats", "mapping_summary.tsv")
INIT_SENTINEL = results_path(".workflow_dirs_initialized")
WORKFLOW_DIRS = [
    results_path("reference"),
    results_path("aligned"),
    results_path("unaligned"),
    results_path("stats"),
    results_path("tmp", "align"),
    log_path("reference"),
    log_path("align"),
    log_path("aligned_index"),
    log_path("stats"),
    log_path("slurm"),
    log_path("workflow"),
    benchmark_path("reference"),
    benchmark_path("align"),
    benchmark_path("aligned_index"),
]


localrules: init_workflow_dirs, summarize_mapping, combine_mapping_stats


rule all:
    input:
        REFERENCE_MMI,
        expand(
            [ALIGNED_BAM, ALIGNED_BAI, UNALIGNED_BAM, MAPPING_STATS],
            sample=SAMPLES,
        ),
        COMBINED_STATS


rule init_workflow_dirs:
    output:
        INIT_SENTINEL
    log:
        log_path("workflow", "init.log")
    conda:
        str(ENV_YAML)
    params:
        dirs=WORKFLOW_DIRS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.dirs:q} "$(dirname {log})"
        touch {output}
        printf "initialized_workflow_dirs\t%s\n" "{output}" > {log}
        """


rule index_rrna_reference:
    input:
        init=INIT_SENTINEL,
        fasta=REFERENCE_FASTA
    output:
        mmi=REFERENCE_MMI
    log:
        log_path("reference", "rrna_35s_index.log")
    benchmark:
        benchmark_path("reference", "rrna_35s_index.tsv")
    threads: config["resources"]["index_threads"]
    conda:
        str(ENV_YAML)
    params:
        minimap2=MINIMAP2_BIN,
        benchmark_dir=str(Path(benchmark_path("reference", "rrna_35s_index.tsv")).parent),
        cluster=rule_cluster_args(
            "index_reference",
            "rrna_mmi",
            log_path("slurm", "rrna_mmi.%A.log"),
            config["resources"]["index_threads"],
        )
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.mmi})" "$(dirname {log})" "{params.benchmark_dir}"
        {params.minimap2} \
          -t {threads} \
          -d {output.mmi} \
          {input.fasta} \
          > {log} 2>&1
        """


rule map_rrna_reads:
    input:
        init=INIT_SENTINEL,
        fastq=sample_fastq,
        index=REFERENCE_MMI
    output:
        aligned=ALIGNED_BAM,
        unaligned=UNALIGNED_BAM
    log:
        log_path("align", "{sample}.log")
    benchmark:
        benchmark_path("align", "{sample}.tsv")
    threads: config["resources"]["align_threads"]
    conda:
        str(ENV_YAML)
    params:
        minimap2=MINIMAP2_BIN,
        preset=MINIMAP2_PRESET,
        extra=MINIMAP2_EXTRA,
        raw_bam=lambda wildcards: results_path("tmp", "align", f"{wildcards.sample}.raw.bam"),
        benchmark_dir=lambda wildcards: str(Path(benchmark_path("align", f"{wildcards.sample}.tsv")).parent),
        cluster=lambda wildcards: rule_cluster_args(
            "align_reads",
            "rrna_map",
            log_path("slurm", f"rrna_map.{wildcards.sample}.%A.log"),
            config["resources"]["align_threads"],
        )
    shell:
        r"""
        set -euo pipefail
        mkdir -p \
          "$(dirname {output.aligned})" \
          "$(dirname {output.unaligned})" \
          "$(dirname {log})" \
          "$(dirname {params.raw_bam})" \
          "{params.benchmark_dir}"
        {params.minimap2} \
          -t {threads} \
          -a \
          -x {params.preset} \
          {params.extra} \
          {input.index} \
          {input.fastq} \
          2> {log} \
        | samtools view -@ {threads} -u - > {params.raw_bam}

        samtools view -@ {threads} -b -F 2308 {params.raw_bam} \
        | samtools sort -@ {threads} -o {output.aligned} -

        samtools view -@ {threads} -b -f 4 -F 2304 {params.raw_bam} \
        | samtools sort -@ {threads} -o {output.unaligned} -

        rm -f {params.raw_bam}
        """


rule index_aligned_bam:
    input:
        init=INIT_SENTINEL,
        bam=ALIGNED_BAM
    output:
        bai=ALIGNED_BAI
    log:
        log_path("aligned_index", "{sample}.log")
    benchmark:
        benchmark_path("aligned_index", "{sample}.tsv")
    threads: config["resources"]["samtools_index_threads"]
    conda:
        str(ENV_YAML)
    params:
        benchmark_dir=lambda wildcards: str(Path(benchmark_path("aligned_index", f"{wildcards.sample}.tsv")).parent),
        cluster=lambda wildcards: rule_cluster_args(
            "index_aligned_bam",
            "bam_index",
            log_path("slurm", f"bam_index.{wildcards.sample}.%A.log"),
            config["resources"]["samtools_index_threads"],
        )
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bai})" "$(dirname {log})" "{params.benchmark_dir}"
        samtools index -@ {threads} {input.bam} {output.bai} > {log} 2>&1
        """


rule summarize_mapping:
    input:
        init=INIT_SENTINEL,
        aligned=ALIGNED_BAM,
        unaligned=UNALIGNED_BAM
    output:
        stats=MAPPING_STATS
    log:
        log_path("stats", "{sample}.log")
    conda:
        str(ENV_YAML)
    script:
        str(SCRIPTS_DIR / "summarize_mapping.py")


rule combine_mapping_stats:
    input:
        expand(MAPPING_STATS, sample=SAMPLES)
    output:
        summary=COMBINED_STATS
    log:
        log_path("stats", "mapping_summary.log")
    conda:
        str(ENV_YAML)
    script:
        str(SCRIPTS_DIR / "combine_mapping_stats.py")
