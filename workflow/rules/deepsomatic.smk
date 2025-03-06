__author__ = "Andrei.G"
__copyright__ = "Copyright 2025, Andrei.G"
__email__ = "andrei.guliaev@scilifelab.uu.se"
__license__ = "GPL-3"


rule deepsomatic_run_deepsomatic:
    input:
        reference=config["reference"]["fasta"],
        normal=os.path.join(config["bam_path"], "{sample}_N.haplotagged.bam"),
        tumor=os.path.join(config["bam_path"], "{sample}_T.haplotagged.bam")
    output:
        vcf="variant_calling/deepsomatic_run_deepsomatic/{sample}.vcf"
    params:
        model=config.get("deepsomatic_run_deepsomatic", {}).get("model", "")
    log:
        dir=directory("variant_calling/deepsomatic_run_deepsomatic/{sample}_run"),
        run="variant_calling/deepsomatic_run_deepsomatic/{sample}.stdout.log"
    benchmark:
        repeat(
            "variant_calling/deepsomatic_run_deepsomatic/{sample}.output.benchmark.tsv",
            config.get("deepsomatic_run_deepsomatic", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("deepsomatic_run_deepsomatic", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("deepsomatic_run_deepsomatic", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepsomatic_run_deepsomatic", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepsomatic_run_deepsomatic", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepsomatic_run_deepsomatic", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepsomatic_run_deepsomatic", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepsomatic_run_deepsomatic", {}).get("container", config["default_container"])
    message:
        "{rule}: do stuff on {wildcards.sample}"
    shell:
        # binding I used during the test run
        # singularity run -B /projects/wp2/nobackup/Twist_Myelom/Test_Analyses/pr_038/analysis/pr_038_001/:/data/
        "run_deepsomatic "
        "--model_type={params.model} "
        "--ref={input.reference} "
        "--reads_normal={input.normal} "
        "--reads_tumor={input.tumor} "
        "--output_vcf={output.vcf} "
        "--sample_name_normal={wildcards.sample}_N "
        "--sample_name_tumor={wildcards.sample}_T "
        "--num_shards={threads} "
        "--logging_dir={log.dir} "
        "--vcf_stats_report=true &> {log.run}"
