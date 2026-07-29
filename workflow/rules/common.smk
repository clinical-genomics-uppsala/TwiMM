__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2025, Andrei Guliaev"
__email__ = "andrei.guliaev@scilifelab.uu.se"
__license__ = "GPL-3"

import itertools
import numpy as np
import pathlib
import pandas as pd
import sys
import yaml
from datetime import datetime
from snakemake.utils import validate
from snakemake.utils import min_version

sys.path.insert(0, str(pathlib.Path(workflow.basedir) / "scripts"))
import utils as pipeline_utils

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

from hydra_genetics.utils.misc import export_config_as_file, get_input_haplotagged_bam, get_input_aligned_bam
from hydra_genetics.utils.software_versions import add_version_files_to_multiqc
from hydra_genetics.utils.software_versions import add_software_version_to_config
from hydra_genetics.utils.software_versions import export_pipeline_version_as_file
from hydra_genetics.utils.software_versions import export_software_version_as_file
from hydra_genetics.utils.software_versions import get_pipeline_version
from hydra_genetics.utils.software_versions import touch_pipeline_version_file_name
from hydra_genetics.utils.software_versions import touch_software_version_file
from hydra_genetics.utils.software_versions import use_container

min_version("7.32.4")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")

try:
    validate(config, schema="../schemas/config.schema.yaml")
except WorkflowError as we:
    # Probably a validation error, but the original exception in lost in
    # snakemake. Pull out the most relevant information instead of a potentially
    # *very* long error message.
    if not we.args[0].lower().startswith("error validating config file"):
        raise
    error_msg = "\n".join(we.args[0].splitlines()[:2])
    parent_rule_ = we.args[0].splitlines()[3].split()[-1]
    if parent_rule_ == "schema:":
        sys.exit(error_msg)
    else:
        schema_hiearachy = parent_rule_.split()[-1]
        schema_section = ".".join(re.findall(r"\['([^']+)'\]", schema_hiearachy)[1::2])
        sys.exit(f"{error_msg} in {schema_section}")

date_string = config.get("pipeline_version", "unknown")
pipeline_version = get_pipeline_version(workflow, pipeline_name="TwiMM")
for pipeline_name in pipeline_version:
    version = pipeline_version[pipeline_name]["version"]
    if version is not None:
        pipeline_version[pipeline_name]["version"] = re.sub(r"[/\\]", "-", version)
version_files = touch_pipeline_version_file_name(pipeline_version, date_string=date_string, directory="results/versions/software")
if use_container(workflow):
    version_files.append(touch_software_version_file(config, date_string=date_string, directory="results/versions/software"))
add_version_files_to_multiqc(config, version_files)


onstart:
    export_pipeline_version_as_file(pipeline_version, date_string=date_string, directory="results/versions/software")
    # Make sure that the user have the requested containers to be used
    if use_container(workflow):
        # From the config retrieve all dockers used and parse labels for software versions. Add
        # this information to config dict.
        update_config, software_info = add_software_version_to_config(config, workflow, False)
        # Print all software used as files. Additional parameters that can be set
        # - directory, default value: software_versions
        # - file_name_ending, default value: mqc_versions.yaml
        # date_string, a string that will be added to the folder name to make it unique (preferably a timestamp)
        export_software_version_as_file(software_info, date_string=date_string, directory="results/versions/software")
        # print config dict as a file. Additional parameters that can be set
        # output_file, default config
        # output_directory, default = None, i.e no folder
        # date_string, a string that will be added to the folder name to make it unique (preferably a timestamp)
    export_config_as_file(update_config, date_string=date_string, directory="results/versions")


### Read and validate resources file

config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")

# Expand {{VAR}} path placeholders in config using top-level string values.
# This allows config.yaml to define e.g. REF_DATA and reference it as {{REF_DATA}}
# in other string values without Snakemake treating REF_DATA as a wildcard.
_cfg_vars = {k: v for k, v in config.items() if isinstance(v, str)}
config.update(pipeline_utils.expand_cfg_placeholders(dict(config), _cfg_vars))
del _cfg_vars

# Format VEP extra config with FASTA path
if "vep" in config and "extra" in config["vep"]:
    config["vep"]["extra"] = config["vep"]["extra"].format(fasta=config["reference"]["fasta"])


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pandas.read_table(config["units"], dtype=str)

if units.platform.iloc[0] in ["PACBIO", "ONT"]:
    units = units.set_index(["sample", "type", "processing_unit", "barcode"], drop=False).sort_index()
else:  # assume that the platform Illumina data with a lane and flowcell columns
    units = units.set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False).sort_index()

validate(units, schema="../schemas/units.schema.yaml")

### Read and validate output file

with open(config["output"]) as output:
    if config["output"].endswith("json"):
        output_spec = json.load(output)
    elif config["output"].endswith("yaml") or config["output"].endswith("yml"):
        output_spec = yaml.safe_load(output.read())

validate(output_spec, schema="../schemas/output_files.schema.yaml")


### Set wildcard constraints
wildcard_constraints:
    sample="|".join(samples.index),
    type="T",


# compile_output_file_list() and generate_copy_rules() are kept here rather than
# moved to workflow/scripts/utils.py: they aren't pure functions of their arguments,
# they read module-level state (output_spec, config, samples, units) and, in
# generate_copy_rules()'s case, register new rules directly on the live `workflow`
# object via workflow.rule/workflow.input/etc. and exec(). That ties them to a real
# Snakemake `workflow` instance at parse time, so they can't be unit tested in
# isolation the way the wildcard-input helper functions below can.
def compile_output_file_list(wildcards):
    outdir = pathlib.Path(output_spec.get("directory", "./"))
    output_files = []

    for f in output_spec["files"]:
        # Please remember to add any additional values down below
        # that the output strings should be formatted with.
        outputpaths = set(
            [
                f["output"].format(sample=sample, type=unit_type)
                for sample in get_samples(samples)
                for unit_type in get_unit_types(units, sample)
            ]
        )

        for op in outputpaths:
            output_files.append(outdir / Path(op))

    return output_files


def generate_copy_rules(output_spec):
    output_directory = pathlib.Path(output_spec.get("directory", "./"))
    rulestrings = []

    for f in output_spec["files"]:
        if f["input"] is None:
            continue

        rule_name = "_copy_{}".format("_".join(re.split(r"\W{1,}", f["name"].strip().lower())))
        input_file = pathlib.Path(f["input"])
        output_file = output_directory / pathlib.Path(f["output"])

        mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
        mem_per_cpu = config.get("_copy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])
        partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
        threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
        time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
        copy_container = config.get("_copy", {}).get("container", config["default_container"])

        rule_code = "\n".join(
            [
                f'@workflow.rule(name="{rule_name}")',
                f'@workflow.input("{input_file}")',
                f'@workflow.output("{output_file}")',
                f'@workflow.log("logs/{rule_name}_{output_file.name}.log")',
                f'@workflow.container("{copy_container}")',
                f'@workflow.resources(time="{time}", threads={threads}, mem_mb="{mem_mb}", '
                f'mem_per_cpu={mem_per_cpu}, partition="{partition}")',
                f'@workflow.shellcmd("{copy_container}")',
                "@workflow.run\n",
                f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, "
                "log, version, rule, conda_env, container_img, singularity_args, use_singularity, "
                "env_modules, bench_record, jobid, is_shell, bench_iteration, cleanup_scripts, "
                "shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
                "__is_snakemake_rule_func=True):",
                '\tshell("(cp --preserve=timestamps {input[0]} {output[0]}) &> {log}", bench_record=bench_record, '
                "bench_iteration=bench_iteration)\n\n",
            ]
        )

        rulestrings.append(rule_code)

    exec(compile("\n".join(rulestrings), "copy_result_files", "exec"), workflow.globals)


generate_copy_rules(output_spec)

# Shallow config copy with aligner replaced by snv_aligner so that
# get_input_aligned_bam returns pbmm2 paths for SNV callers while SV callers
# and whatshap_haplotag continue to use the primary (vacmap) aligner.
_snv_config = {**config, "aligner": config.get("snv_aligner", config["aligner"])}


# The functions below are thin wrappers around workflow/scripts/utils.py: the
# actual logic lives there (as plain functions of explicit arguments) so it can be
# unit tested without a live Snakemake `workflow`/`config`/`units` context. These
# wrappers just supply the module-level config/units that Snakemake's `input:`/
# `params:` callables (which only receive `wildcards`) can't be handed directly.
def get_local_vcfs_for_svdb_merge(wildcards, add_suffix=False):
    return pipeline_utils.get_local_vcfs_for_svdb_merge(wildcards, config, add_suffix=add_suffix)


def get_svdb_merge_priority(wildcards):
    return pipeline_utils.get_svdb_merge_priority(wildcards, config)


def get_tc_file(wildcards):
    return pipeline_utils.get_tc_file(wildcards, config)


def get_snv_caller_output(wildcards):
    return pipeline_utils.get_snv_caller_output(wildcards, config)


get_concat_caller_vcfs = pipeline_utils.get_concat_caller_vcfs


def get_ubam_input(wildcards):
    return pipeline_utils.get_ubam_input(wildcards, units)
