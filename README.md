# Strelka2 VCF Annotator

A command-line Python tool to annotate somatic (tumor/normal) SNV and indel VCFs produced by Strelka2 with additional metrics like ADs, VAFs, and Fisher test statistics.
It supports merged VCFs, and VEP-annotated outputs.


## ✨ Features

- Works with one or more VCF files (SNVs and/or indels)
- Merges VCFs on-the-fly if multiple are given
- Supports filtering to PASS variants only using `bcftools`
- Adds both simple and total read depth-based VAFs (discussed below)
- Computes Fisher exact test statistics
- Optionally indexes the output with `bcftools`


## 🔧 Installation and dependencies

- Python 3
- Python packages:
  - `scipy`
- `bcftools`

We recommend using a conda environment for reproducibility.

```bash
# Create a new environment with Python
conda create -n strelka2_vcf_annotator python=3.10 -y
conda activate strelka2_vcf_annotator

# Install bcftools
conda install -c bioconda bcftools

# Install scipy
conda install scipy
```

Alternatively, install `bcftools` manually:
- [Build and install from source](https://www.htslib.org/download/)

Ensure `bcftools` is in your `$PATH`:
```bash
which bcftools
```


## 🚀 Usage

### 🔄 Clone the repository

First, clone this repository and enter the directory:
```bash
git clone https://github.com/blazoder/strelka2-vcf-annotator.git
cd strelka2-vcf-annotator
```

Then run the script as shown below:

### Simple:
```bash
python annotate_strelka2.py \
  --input sample1.vcf.gz \
  --output sample1.annotated.vcf.gz
```

### With merging:
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


## 🧬 Fields added

### INFO fields added:
- `TUMREF`, `TUMALT`, `NORMREF`, `NORMALT`
- `TUMVAF`: ALT / (REF + ALT)
- `TUMVAF_TOTAL`: ALT / total depth (all bases)
- `TUMVARFRACTION`: ALT in tumor / total ALT (in tumor + normal)
- `LOG_FISHER`: Fisher's test on REF vs ALT
- `LOG_FISHER_TOTAL`: Fisher's test on ALT vs all bases

### FORMAT fields added:
- `AD`: REF and ALT counts (comma-separated)
- `DPVAF`: REF + ALT
- `DPVAF_TOTAL`: Total depth (all bases)
- `VAF`: ALT / (REF + ALT)
- `VAF_TOTAL`: ALT / total depth (all bases)

> Tier1 counts only are considered. Read more about Strelka2 [here](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#vcf-files).

`TUMVAF` and `LOG_FISHER` from INFO, and `DPVAF`, `VAF` from FORMAT are calculated using only the REF and ALT alleles, in line with Strelka2's recommended approach.

However, sequencing reads may also support other nucleotides or alleles at the same position. To provide a more comprehensive view of allelic representation, we additionally compute:

- `TUMVAF_TOTAL` and `LOG_FISHER_TOTAL` (in INFO)
- `VAF_TOTAL` and `DPVAF_TOTAL` (in FORMAT)

These fields use the total tier1 depth from all observed bases, not just the REF and ALT alleles. This helps better reflect situations with mixed or noisy signal beyond the canonical REF/ALT model.


## 🗂 Output
The final output is a compressed VCF (`.vcf.gz`) with all annotations added. If `--index-output` is used, it also generates a `.tbi` index file.


## 📝 License
This tool is released under the [MIT License](LICENSE). You are free to use, modify, and distribute it with proper attribution.

---
MIT License © 2025
Blaz Oder | Clinical Genetics | Karolinska Institutet
