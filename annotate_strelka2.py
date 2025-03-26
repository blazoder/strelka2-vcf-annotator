#!/usr/bin/env python3

import sys
import argparse
import gzip
import math
import subprocess
import tempfile
from pathlib import Path
from scipy.stats import fisher_exact

def calculate_log_fisher(tumor_ref, tumor_alt, normal_ref, normal_alt):
    table = [[tumor_ref, tumor_alt], [normal_ref, normal_alt]]
    _, p = fisher_exact(table)
    return -math.log10(p) if p > 0 else 300.0

def calculate_log_fisher_total(tumor_alt, tumor_total, normal_alt, normal_total):
    tumor_other = tumor_total - tumor_alt
    normal_other = normal_total - normal_alt
    table = [[tumor_other, tumor_alt], [normal_other, normal_alt]]
    _, p = fisher_exact(table)
    return -math.log10(p) if p > 0 else 300.0

def calculate_tumvarfraction(tum_alt, norm_alt):
    total_alt = tum_alt + norm_alt
    return tum_alt / total_alt if total_alt > 0 else 0.0

def parse_counts(format_keys, sample_values, ref_base=None, alt_base=None, variant_type='snv'):
    data = dict(zip(format_keys, sample_values.split(":")))
    if variant_type == 'snv':
        ref_key = ref_base + "U"
        alt_key = alt_base + "U"
        ref_counts = data.get(ref_key, "0,0").split(",")
        alt_counts = data.get(alt_key, "0,0").split(",")
        total_tier1 = sum(int(data.get(b + "U", "0,0").split(",")[0]) for b in "ACGT" if b + "U" in data)
    else:
        ref_counts = data.get("TAR", "0,0").split(",")
        alt_counts = data.get("TIR", "0,0").split(",")
        total_tier1 = int(ref_counts[0]) + int(alt_counts[0])

    return int(ref_counts[0]), int(alt_counts[0]), total_tier1

def is_snv(ref, alt):
    return len(ref) == 1 and len(alt) == 1

def open_vcf(file_path):
    return gzip.open(file_path, "rt") if file_path.endswith(".gz") else open(file_path, "r")

def write_vcf(output_path):
    return gzip.open(output_path, "wt") if output_path.endswith(".gz") else open(output_path, "w")

def filter_pass_variants(input_vcf):
    filtered_vcf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf.gz").name
    subprocess.run(["bcftools", "view", "-f", "PASS", "-O", "z", "-o", filtered_vcf, input_vcf], check=True)
    subprocess.run(["bcftools", "index", "-f", filtered_vcf], check=True)
    return filtered_vcf

def annotate_vcf(input_vcf, output_vcf, only_pass=False, index_output=False):
    if only_pass:
        input_vcf = filter_pass_variants(input_vcf)

    with open_vcf(input_vcf) as fin, write_vcf(output_vcf) as fout:
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
            elif line.startswith("#CHROM"):
                fout.write('##INFO=<ID=TUMREF,Number=1,Type=Integer,Description="Tumor REF read count (tier1)">\n')
                fout.write('##INFO=<ID=TUMALT,Number=1,Type=Integer,Description="Tumor ALT read count (tier1)">\n')
                fout.write('##INFO=<ID=NORMREF,Number=1,Type=Integer,Description="Normal REF read count (tier1)">\n')
                fout.write('##INFO=<ID=NORMALT,Number=1,Type=Integer,Description="Normal ALT read count (tier1)">\n')
                fout.write('##INFO=<ID=TUMVAF,Number=1,Type=Float,Description="Tumor VAF = ALT / (REF + ALT)">\n')
                fout.write('##INFO=<ID=TUMVAF_TOTAL,Number=1,Type=Float,Description="Tumor VAF = ALT / total tier1 depth">\n')
                fout.write('##INFO=<ID=TUMVARFRACTION,Number=1,Type=Float,Description="ALT reads in tumor / total ALT reads">\n')
                fout.write('##INFO=<ID=LOG_FISHER,Number=1,Type=Float,Description="-log10(p) Fisher test on REF/ALT">\n')
                fout.write('##INFO=<ID=LOG_FISHER_TOTAL,Number=1,Type=Float,Description="-log10(p) Fisher test on ALT vs rest (tier1)">\n')
                fout.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for REF,ALT (tier1)">\n')
                fout.write('##FORMAT=<ID=DPVAF,Number=1,Type=Integer,Description="Depth = REF + ALT">\n')
                fout.write('##FORMAT=<ID=DPVAF_TOTAL,Number=1,Type=Integer,Description="Total tier1 depth (all bases)">\n')
                fout.write('##FORMAT=<ID=VAF,Number=A,Type=Float,Description="ALT / (REF + ALT)">\n')
                fout.write('##FORMAT=<ID=VAF_TOTAL,Number=A,Type=Float,Description="ALT / total tier1 depth">\n')
                fout.write(line)
            else:
                parts = line.strip().split("\t")
                ref, alt = parts[3], parts[4]
                fmt_keys = parts[8].split(":")
                normal_sample = parts[9]
                tumor_sample = parts[10]

                vtype = 'snv' if is_snv(ref, alt) else 'indel'

                norm_ref, norm_alt, norm_total = parse_counts(fmt_keys, normal_sample, ref, alt, vtype)
                tum_ref, tum_alt, tum_total = parse_counts(fmt_keys, tumor_sample, ref, alt, vtype)

                tumvaf = tum_alt / (tum_alt + tum_ref) if (tum_alt + tum_ref) > 0 else 0.0
                normvaf = norm_alt / (norm_alt + norm_ref) if (norm_alt + norm_ref) > 0 else 0.0
                tumvaf_total = tum_alt / tum_total if tum_total > 0 else 0.0
                normvaf_total = norm_alt / norm_total if norm_total > 0 else 0.0

                tumvarfraction = calculate_tumvarfraction(tum_alt, norm_alt)
                log_fisher = calculate_log_fisher(tum_ref, tum_alt, norm_ref, norm_alt)
                log_fisher_total = calculate_log_fisher_total(tum_alt, tum_total, norm_alt, norm_total)

                info = parts[7]
                extra_info = f"TUMREF={tum_ref};TUMALT={tum_alt};NORMREF={norm_ref};NORMALT={norm_alt};"
                extra_info += f"TUMVAF={tumvaf:.4f};TUMVAF_TOTAL={tumvaf_total:.4f};"
                extra_info += f"TUMVARFRACTION={tumvarfraction:.4f};LOG_FISHER={log_fisher:.4f};LOG_FISHER_TOTAL={log_fisher_total:.4f}"
                parts[7] = info + ";" + extra_info

                for tag in ["AD", "DPVAF", "DPVAF_TOTAL", "VAF", "VAF_TOTAL"]:
                    if tag not in fmt_keys:
                        fmt_keys.append(tag)
                parts[8] = ":".join(fmt_keys)

                normal_fields = normal_sample.split(":") + [f"{norm_ref},{norm_alt}", str(norm_ref + norm_alt), str(norm_total), f"{normvaf:.4f}", f"{normvaf_total:.4f}"]
                tumor_fields = tumor_sample.split(":") + [f"{tum_ref},{tum_alt}", str(tum_ref + tum_alt), str(tum_total), f"{tumvaf:.4f}", f"{tumvaf_total:.4f}"]
                parts[9] = ":".join(normal_fields)
                parts[10] = ":".join(tumor_fields)

                fout.write("\t".join(parts) + "\n")

    if index_output:
        subprocess.run(["bcftools", "index", "-f", output_vcf])

def merge_vcfs(vcf_list):
    temp_output = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf.gz").name
    subprocess.run(["bcftools", "concat", "-a", "-O", "z", "-o", temp_output] + vcf_list, check=True)
    return temp_output

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate Strelka2 VCF(s) with AD, VAF, total VAF, and Fisher scores")
    parser.add_argument("--input", nargs='+', required=True, help="One or more input VCFs (SNVs/INDELs)")
    parser.add_argument("--output", required=True, help="Output annotated VCF")
    parser.add_argument("--only-pass", action="store_true", help="Filter variants to PASS only using bcftools")
    parser.add_argument("--index-output", action="store_true", help="Index final output with bcftools")
    args = parser.parse_args()

    if len(args.input) == 1:
        input_vcf = args.input[0]
    else:
        input_vcf = merge_vcfs(args.input)

    annotate_vcf(input_vcf, args.output, only_pass=args.only_pass, index_output=args.index_output)
