__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2025, Andrei Guliaev"
__email__ = "andrei.guliae@scilifelab.uu.se"
__license__ = "GPL-3"


rule pbsv_call_structural_variants:
    input:
        bam=os.path.join(config["bam_path"], "{sample}_T.haplotagged.bam"),
    output:
        svsig="variant_calling/pbsv_call_structural_variants/{sample}.svsig.gz",
        vcf="variant_calling/pbsv_call_structural_variants/{sample}.vcf"
    params:
        trf=config.get("pbsv_call_structural_variants", {}).get("trf", ""),
        reference=config.get("pbsv_call_structural_variants", {}).get("reference", ""),
        extra=config.get("pbsv_call_structural_variants", {}).get("extra", ""),
    log:
        discover="variant_calling/pbsv_call_structural_variants/{sample}.discover.log",
        call="variant_calling/pbsv_call_structural_variants/{sample}.call.log"
    benchmark:
        repeat(
            "variant_calling/pbsv_call_structural_variants/{sample}.output.benchmark.tsv",
            config.get("pbsv_call_structural_variants", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("pbsv_call_structural_variants", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("pbsv_call_structural_variants", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pbsv_call_structural_variants", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pbsv_call_structural_variants", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("pbsv_call_structural_variants", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pbsv_call_structural_variants", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("pbsv_call_structural_variants", {}).get("container", config["default_container"])
    message:
        "{rule}: do stuff on {input.bam}"
    shell:
        "pbsv discover "
        "--tandem-repeats {params.trf} "
        "--log-file {log.discover} "
        "{input.bam} {output.svsig} && "
        "pbsv call -j {threads} "
        "--log-file {log.call} "
        "--hifi {params.reference} "
        "{output.svsig} "
        "{output.vcf}"
