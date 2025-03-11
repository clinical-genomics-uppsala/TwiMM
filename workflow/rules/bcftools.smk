__author__ = "andrei guliaev"
__copyright__ = "Copyright 2025, andrei guliaev"
__email__ = "andrei.guliaev@scilifielab.uu.se"
__license__ = "GPL-3"


rule bcftools_filter_include_region:
    input:
        vcf="variant_calling/deepsomatic_run_deepsomatic/{sample}.vcf",
    output:
        vcf="variant_calling/filtered/{sample}.include.panel.vcf",
    params:
        bed=config["bcftools_filter_include_region"]["panel"]
    log:
        "filtering/bcftools_filter_include_region/{sample}.output.log",
    benchmark:
        repeat(
            "filtering/bcftools_filter_include_region/{sample}.output.benchmark.tsv",
            config.get("bcftools_filter_include_region", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("bcftools_filter_include_region", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcftools_filter_include_region", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_filter_include_region", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_filter_include_region", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_filter_include_region", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_filter_include_region", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_filter_include_region", {}).get("container", config["default_container"])
    message:
        "{rule}: do stuff on {input.vcf}"
    shell:
        "(bcftools filter "
        "-R {params.bed} "
        "{input.vcf} "
        "-o {output.vcf}) &> {log}"
