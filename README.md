🧬 strelka-vcf-annotator

A lightweight Python tool to annotate Strelka2 somatic VCFs with tumor/normal read counts, variant allele frequencies (VAF), and Fisher’s exact test p-values.

If you provide multiple input VCFs, they will be automatically merged before annotation using bcftools.

🚀 Features

✅ Annotates tumor and normal:
REF and ALT read counts
Variant Allele Frequencies (VAF)
Fisher’s exact test p-values (LOG_FISHER)
Tumor ALT fraction (TUMVARFRACTION)
✅ Adds new FORMAT fields:
AD — Allelic depths (REF, ALT)
DPVAF — Depth used for VAF
VAF — Variant Allele Frequency
✅ Supports:
One or more input VCFs (SNVs and/or INDELs)
Automatic merging when multiple files are given
Filtering to PASS variants only (--only-pass)
Optional indexing of the final output VCF (--index-output)
VEP-annotated VCFs (annotations are preserved)
📦 Installation

Requires Python 3 and scipy.

Install Python dependencies:

pip install scipy
Install bcftools (required for merging and indexing):

conda install -c bioconda bcftools
🛠️ Usage

1. Annotate a single VCF file
python3 strelka_vcf_annotator.py \
  --input somatic.vcf \
  --output somatic.annotated.vcf
2. Annotate and keep only PASS variants
python3 strelka_vcf_annotator.py \
  --input somatic.vcf \
  --output somatic.annotated.vcf \
  --only-pass
3. Annotate and automatically merge SNV + INDEL VCFs
python3 strelka_vcf_annotator.py \
  --input snvs.vcf indels.vcf \
  --output merged.annotated.vcf
4. Annotate, filter to PASS, and index the result
python3 strelka_vcf_annotator.py \
  --input snvs.vcf indels.vcf \
  --output merged.annotated.vcf \
  --only-pass \
  --index-output
Note: If more than one file is given to --input, the script will merge them automatically using bcftools.

🧪 Output Annotations

INFO fields added
Field	Description
TUMREF	Tumor REF count (tier1)
TUMALT	Tumor ALT count (tier1)
NORMREF	Normal REF count (tier1)
NORMALT	Normal ALT count (tier1)
TUMVAF	Tumor ALT / (REF + ALT)
TUMVARFRACTION	Tumor ALT / (Tumor ALT + Normal ALT)
LOG_FISHER	-log10(p-value) from Fisher’s exact test
FORMAT fields added
Field	Description
AD	Allelic depths: REF,ALT
DPVAF	Depth at site used to calculate VAF
VAF	ALT / (REF + ALT)
💬 Example

Before annotation:

#CHROM  POS   ID  REF ALT QUAL FILTER INFO FORMAT     NORMAL        TUMOR
chr1    12345 .   A   T   .    PASS   ...  AU:CU:GU    ...           ...
After annotation:

#CHROM  POS   ID  REF ALT QUAL FILTER INFO                                                                 FORMAT                             NORMAL                     TUMOR
chr1    12345 .   A   T   .    PASS  TUMREF=30;TUMALT=15;TUMVAF=0.333;...                                  AU:CU:GU:AD:DPVAF:VAF              ...:30,2:32:0.063           ...:35,15:50:0.300
📄 License

MIT License © 2024
Blaz Oder | Clinical Genetics | Karolinska Institutet
# Strelka2 VCF Annotator

A command-line Python tool to annotate SNV and indel VCFs produced by Strelka2 with additional metrics like VAFs and Fisher test statistics.

---

## ✨ Features

- Works with one or more VCF files (SNVs and/or indels)
- Merges VCFs on-the-fly if multiple are given
- Supports filtering to PASS variants only using `bcftools`
- Adds both simple and total read depth-based VAFs
- Computes Fisher exact test statistics
- Optionally indexes the output with `bcftools`

---

## 🔧 Dependencies

- Python 3
- `bcftools` (e.g. via conda or modules)
- Python packages:
  - `scipy`

Install dependencies with:
```bash
pip install scipy
# and ensure bcftools is in your $PATH
```

---

## 🚀 Usage

```bash
python annotate_strelka2.py \
  --input sample1_snv.vcf.gz sample1_indel.vcf.gz \
  --output sample1.annotated.vcf.gz
```

### With PASS filtering:
```bash
python annotate_strelka2.py \
  --input sample1_snv.vcf.gz sample1_indel.vcf.gz \
  --output sample1.pass_only.annotated.vcf.gz \
  --only-pass
```

### With indexing of output:
```bash
python annotate_strelka2.py \
  --input sample1_snv.vcf.gz sample1_indel.vcf.gz \
  --output sample1.annotated.vcf.gz \
  --only-pass \
  --index-output
```

---

## 🧬 Fields Annotated

### INFO fields added:
- `TUMREF`, `TUMALT`, `NORMREF`, `NORMALT`: Tier1 counts
- `TUMVAF`: ALT / (REF + ALT)
- `TUMVAF_TOTAL`: ALT / total tier1 depth
- `TUMVARFRACTION`: ALT in tumor / total ALT
- `LOG_FISHER`: Fisher's test on REF vs ALT
- `LOG_FISHER_TOTAL`: Fisher's test on ALT vs all other reads

### FORMAT fields added:
- `AD`: REF and ALT tier1 counts
- `DPVAF`: REF + ALT depth
- `DPVAF_TOTAL`: Total tier1 depth
- `VAF`: ALT / (REF + ALT)
- `VAF_TOTAL`: ALT / total tier1 depth

---

## 🗂 Output
The final output is a compressed VCF (`.vcf.gz`) with all annotations added. If `--index-output` is used, it also generates a `.tbi` index file.

---

## 📫 Contact
Created by Blaz for high-throughput post-Strelka2 variant processing. Reach out for suggestions, issues or feature requests!

