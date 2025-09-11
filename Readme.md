# Variant Calling and Annotation

This repository demonstrates a workflow for **variant discovery** and **functional annotation** using widely adopted bioinformatics tools.

* **Variant Calling:** Performed with [GATK (Genome Analysis Toolkit)](https://gatk.broadinstitute.org/).
* **Annotation:** Performed with [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/).

The pipeline starts with aligned sequencing data and outputs annotated variants ready for downstream analysis.

---

## 📂 Repository Contents

* `data/` – Sample input data for testing the pipeline
* `scripts/` – Shell/Python scripts for running GATK and ANNOVAR
* `results/` – Example outputs from the workflow
* `README.md` – Project overview and usage instructions

---

## 🚀 Workflow Overview

1. **Preprocessing**

  * QC: Fastqc
   * Trimming and adapter removal
   * Input: aligned sequencing data (BAM format)
   * Steps: marking duplicates, base quality score recalibration (BQSR)

1. **Variant Calling (GATK)**

   * Tools: `HaplotypeCaller`, `GenotypeGVCFs`
   * Output: raw VCF file

2. **Annotation (ANNOVAR)**

   * Input: VCF from GATK
   * Steps: functional annotation, gene-based annotation, and database-driven annotation (dbSNP, ClinVar, etc.)
   * Output: annotated variant tables

---

## ⚙️ Requirements

* **Software:**

  * [GATK ≥ 4.x](https://gatk.broadinstitute.org/)
  * [ANNOVAR](https://annovar.openbioinformatics.org/)
  * [samtools](http://www.htslib.org/)


* **Data:**

  * Reference genome (FASTA)
  * Aligned reads (BAM)
  * VCF
  * Annotated VCF

## 📜 License

This project is released for **educational and research use**. Please check the license terms of **GATK** and **ANNOVAR** before commercial use.

---

