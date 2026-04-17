__author__ = "Martin Rippin"
__copyright__ = "Copyright 2021, Martin Rippin"
__email__ = "martin.rippin@igp.uu.se"
__license__ = "GPL-3"

# Local copy of hydra-genetics/snv_indels fix_af.py with snv_concat added as a
# valid caller (treated identically to clairs_to: FORMAT/AF is propagated to
# INFO/AF with no transformation). All other logic is unchanged from upstream.

import logging

import pysam


def fixFreebayes(header, row):
    ad = row.samples[0].get("AD")
    if ad is None:
        return None
    total = sum(ad)
    if total == 0:
        return None
    return tuple(a / total for a in ad[1:])


def getCaller(path: str):
    pathParts = path.split("/")
    if len(pathParts) == 3:
        return pathParts[1]
    else:
        raise ValueError(
            "{} is not a valid input for this script. Required:"
            "'snv_indels/caller/sample_type.vcf'.".format(path)
        )


def modifyHeader(caller: str, header: pysam.libcbcf.VariantHeader):
    if caller in ("pisces", "freebayes", "pbrun_deepvariant", "deepvariant", "varscan"):
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", "AF"),
                ("Number", "A"),
                ("Type", "Float"),
                ("Description", "Allele frequency"),
            ],
        )
    if caller in (
        "pisces",
        "pbrun_mutectcaller_T",
        "gatk_mutect2",
        "pbrun_deepvariant",
        "deepvariant",
        "clairs_to",
        "snv_concat",  # concat of clairs_to + deepsomatic; FORMAT/AF already present
    ):
        header.info.add("AF", "A", "Float", "Allele frequency")
    elif caller == "varscan":
        header.info.add(
            "AF",
            "A",
            "Float",
            "Allele count divided on depth (Quality of bases: Phred score >= 15)",
        )
    return header


def writeNewVcf(
    path: str,
    header: pysam.libcbcf.VariantHeader,
    vcf: pysam.libcbcf.VariantFile,
    caller: str,
):
    new_vcf = pysam.VariantFile(path, "w", header=header)
    for row in vcf.fetch():
        if caller == "freebayes":
            row.info["AF"] = fixFreebayes(header, row)
            row.samples[0]["AF"] = fixFreebayes(header, row)
        elif caller == "haplotypecaller":
            row.info["AF"] = row.samples[0].get("AF")
        elif caller == "gatk_mutect2":
            row.info["AF"] = row.samples[0].get("AF")
        elif caller == "pisces":
            row.info["AF"] = row.samples[0].get("VF")
            row.samples[0]["AF"] = row.samples[0].get("VF")
        elif caller == "vardict":
            row.info["AF"] = row.samples[0].get("AF")
        elif caller == "pbrun_mutectcaller_T":
            row.info["AF"] = row.samples[0].get("AF")
        elif caller == "pbrun_deepvariant":
            row.info["AF"] = row.samples[0].get("VAF")
            row.samples[0]["AF"] = row.samples[0].get("VAF")
        elif caller == "deepvariant":
            row.info["AF"] = row.samples[0].get("VAF")
            row.samples[0]["AF"] = row.samples[0].get("VAF")
        elif caller in ("clairs_to", "snv_concat"):
            # ClairS-TO reports FORMAT/AF; DeepSomatic reports FORMAT/VAF.
            # Try AF first, fall back to VAF so both caller types are handled.
            af = row.samples[0].get("AF")
            if af is None:
                af = row.samples[0].get("VAF")
            if af is not None:
                row.info["AF"] = af
            else:
                logging.warning(
                    "Record at %s:%s has neither FORMAT/AF nor FORMAT/VAF; INFO/AF will not be set.",
                    row.chrom,
                    row.pos,
                )
        elif caller == "gatk_select_variants_final":
            row.info["AF"] = row.samples[0].get("AF")
        elif caller == "varscan":
            row.info["AF"] = row.samples[0].get("AD") / row.samples[0].get("DP")
            row.samples[0]["AF"] = row.samples[0].get("AD") / row.samples[0].get("DP")
        else:
            raise ValueError(
                "{} is not a valid caller for this script. Choose between: "
                "freebayes, haplotypecaller, gatk_mutect2, gatk_select_variants_final, "
                "pbrun_deepvariant, deepvariant, clairs_to, pisces, vardict, varscan, "
                "snv_concat.".format(caller)
            )
        new_vcf.write(row)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, filename=snakemake.log[0])  # noqa: F821
    logging.info("Read %s", snakemake.input.vcf)  # noqa: F821
    vcf = pysam.VariantFile(snakemake.input.vcf)  # noqa: F821
    logging.info("Determine caller...")
    caller = getCaller(snakemake.input.vcf)  # noqa: F821
    logging.info("Caller is %s", caller)
    logging.info("Add info to header if necessary")
    header = modifyHeader(caller, vcf.header)
    logging.info("Start writing to %s", snakemake.output.vcf)  # noqa: F821
    writeNewVcf(snakemake.output.vcf, caller=caller, header=header, vcf=vcf)  # noqa: F821
    logging.info("Successfully written vcf file %s", snakemake.output.vcf)  # noqa: F821
