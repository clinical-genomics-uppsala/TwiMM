__author__ = "Andrei Guliaev"
__copyright__ = "Copyright 2025, Andrei Guliaev"
__email__ = "andrei.guliaev@scilifelab.uu.se"
__license__ = "GPL-3"

import re


def expand_cfg_placeholders(obj, cfg_vars):
    """Expand {{VAR}} placeholders in obj using top-level string config values.

    This allows config.yaml to define e.g. REF_DATA and reference it as {{REF_DATA}}
    in other string values without Snakemake treating REF_DATA as a wildcard.
    """
    if isinstance(obj, str):

        def _replace(m):
            key = m.group(1)
            if key not in cfg_vars:
                raise KeyError(f"Config placeholder '{{{{{key}}}}}' has no matching top-level string variable.")
            return cfg_vars[key]

        return re.sub(r"\{\{(\w+)\}\}", _replace, obj)
    if isinstance(obj, list):
        return [expand_cfg_placeholders(i, cfg_vars) for i in obj]
    if isinstance(obj, dict):
        return {k: expand_cfg_placeholders(v, cfg_vars) for k, v in obj.items()}
    return obj


def get_local_vcfs_for_svdb_merge(wildcards, config, add_suffix=False):
    """Return native SV caller VCF paths for SVDB merge for the given tc_method wildcard.
    When add_suffix=False (default): returns plain file paths suitable for
    use as Snakemake input entries.
    When add_suffix=True: appends ':{caller}' to each path, producing the
    'file.vcf:CALLER' notation expected by SVDB --vcf; do NOT use these
    strings as Snakemake input files.
    """
    if "svdb_merge" not in config:
        raise KeyError("Config section 'svdb_merge' is missing.")
    if "tc_method" not in config["svdb_merge"]:
        raise KeyError("Config section 'svdb_merge: tc_method' is missing.")

    for v in config["svdb_merge"]["tc_method"]:
        if "name" not in v:
            raise KeyError("A 'tc_method' entry in config['svdb_merge'] is missing the 'name' key.")
        if v["name"] != wildcards.tc_method:
            continue

        if "sv_caller" not in v:
            raise KeyError(f"tc_method '{wildcards.tc_method}' is missing the 'sv_caller' key.")

        vcf_paths = []
        for entry in v["sv_caller"]:
            if "caller" not in entry or "vcf" not in entry:
                raise KeyError(f"Each sv_caller entry must have 'caller' and 'vcf' keys, got: {entry}")
            caller_suffix = f":{entry['caller']}" if add_suffix else ""
            vcf_paths.append(entry["vcf"].format(sample=wildcards.sample, type=wildcards.type) + caller_suffix)
        return vcf_paths

    available = [v["name"] for v in config["svdb_merge"]["tc_method"] if "name" in v]
    raise KeyError(f"tc_method '{wildcards.tc_method}' not found in 'svdb_merge' config. Available: {available}")


def get_svdb_merge_priority(wildcards, config):
    if "svdb_merge" not in config:
        raise KeyError("Config section 'svdb_merge' is missing.")
    if "tc_method" not in config["svdb_merge"]:
        raise KeyError("Config section 'svdb_merge: tc_method' is missing.")

    priority = None
    for v in config["svdb_merge"]["tc_method"]:
        if "name" not in v:
            raise KeyError("A 'tc_method' entry in config['svdb_merge'] is missing the 'name' key.")
        tc_method_name = v["name"]

        if "sv_caller" not in v:
            raise KeyError(f"tc_method '{tc_method_name}' is missing the 'sv_caller' key.")

        if "priority" not in v:
            raise KeyError(f"tc_method '{tc_method_name}' is missing the 'priority' key.")

        if tc_method_name == wildcards.tc_method:
            priority = v["priority"]

    if priority is not None:
        return priority

    raise KeyError(f"tc_method '{wildcards.tc_method}' not found in 'svdb_merge' configuration.")


def get_tc_file(wildcards, config):
    tc_method = wildcards.tc_method
    if tc_method == "pathology":
        return config["samples"]
    else:
        return f"cnv_sv/{tc_method}_purity_file/{wildcards.sample}_{wildcards.type}.purity.txt"


def get_snv_caller_output(wildcards, config):
    if config["use_deepsomatic"]:
        # Both callers are used; fix_af receives the sorted concat VCF
        return "snv_indels/snv_concat/{sample}_{type}.vcf.gz"
    else:
        return "snv_indels/clairs_to/{sample}_{type}.vcf.gz"


def get_concat_caller_vcfs(wildcards):
    """Return caller-tagged VCF paths for bcftools concat.
    Only called when use_deepsomatic is true (the concat rule is only
    triggered by get_snv_caller_output in that case).
    """
    return [
        f"snv_indels/clairs_to/{wildcards.sample}_{wildcards.type}.caller_tagged.vcf.gz",
        f"snv_indels/deepsomatic_t_only/{wildcards.sample}_{wildcards.type}.caller_tagged.vcf.gz",
    ]


def get_ubam_input(wildcards, units):
    return units.loc[(wildcards.sample, wildcards.type), "bam"].tolist()
