# Pipeline setup and configuration
There are a number of main files that rule the pipeline's execution:

* Snakefile
* common.smk
* config.yaml
* resources.yaml
* samples.tsv
* units.tsv

There is more general information about the content of these files in hydra-genetics documentation in [code standards](https://hydra-genetics.readthedocs.io/en/latest/development/standards/), [config](https://hydra-genetics.readthedocs.io/en/latest/make_pipeline/config/) and [Snakefile](https://hydra-genetics.readthedocs.io/en/latest/make_pipeline/import/).

## Snakefile
The `Snakefile` is located in workflow/ and imports hydra-genetics modules and rules as well as modifies these rules when needed. 
It also imports pipeline specific rules and define rule orders. Finally, this is where the `rule all` is defined.

## common.smk
The `common.smk` is located under workflow/rules/. 
This is a general rule taking care of any actions that are not directly connected with running a specific program. 
It includes version checks, import of config, resources, tsv-files and validations using schemas. 
Functions used by pipeline specific rules are also defined here as well as the output files using the function **compile_output_list** which programmatically generates a list of all necessary output files for the module to be targeted in the all rule defined in the `Snakemake` file. 
See further [Result files](https://hydra-genetics.readthedocs.io/en/latest/make_pipeline/results/).

## config.yaml
The `config.yaml` is located under config/. 
The file ties all the inputs, outputs and other dependencies as well as parameters for different rules together.
See further [pipeline configuration](https://hydra-genetics.readthedocs.io/en/latest/make_pipeline/config/).

## resources.yaml
The `resources.yaml` is located under config/. 
The file declares default resources used by rules as well as resources for specific rules that need more than allocated by default. 
See further [pipeline configuration](https://hydra-genetics.readthedocs.io/en/latest/make_pipeline/config/).

```yaml
# ex, default resources
default_resources:
  threads: 1
  time: "4:00:00"
  mem_mb: 6144
  mem_per_cpu: 6144
  partition: "low"
```


## profile.yaml
Profiles are saved as yaml files to control how snakemake will be executed, if jobs will be submitted
to a cluster, use singularity, restart on failure and so forth. It also forwards requested resources to *drmaa* using
a drmaa variable.

```yaml
# ex, snakemake settings
jobs: 100
keep-going: True
restart-times: 2
rerun-incomplete: True
use-singularity: True
configfile: "config/config.yaml"
singularity-args: "-e --cleanenv -B /projects -B /data -B /beegfs
```

```yaml
# ex, drmaa settings
drmaa: " -A wp1 -N 1-1 -t {resources.time} -n {resources.threads} --mem={resources.mem_mb} --mem-per-cpu={resources.mem_per_cpu} --mem-per-cpu={resources.mem_per_cpu} --partition={resources.partition} -J {rule} -e slurm_out/{rule}_%j.err -o slurm_out/{rule}_%j.out"
drmaa-log-dir: "slurm_out"
default-resources: [threads=1, time="04:00:00", partition="low", mem_mb="3074", mem_per_cpu="3074"]
```

## samples.tsv and units.tsv
The `samples.tsv` and `units.tsv` are input files that must be generated before running the pipeline and should in general be located in the base folder of the analysis folder, can be changed in the config.yaml. See further [running the pipeline](running.md) and [create input files](https://hydra-genetics.readthedocs.io/en/latest/run_pipeline/create_sample_files/).

### Example of `samples.tsv`

| sample | tumor_content |
|--------|---------------|
| MM0001 | 0.8           |
| MM0002 | 0.8           |


### Example of `units.tsv`

| sample | type | platform | machine | processing_unit         | barcode | methylation | bam                 |
|--------|------|----------|---------|-------------------------|---------|-------------|---------------------|
| MM0001 | T    | PACBIO   | REVIO   | m84045_031032_151755_s4 | XXX     | No          | /path/to/MM0001.bam |
| MM0002 | T    | PACBIO   | REVIO   | m84045_031032_151756_s5 | ZZZ     | IDK         | /path/to/MM0002.bam |