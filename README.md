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

