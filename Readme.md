# Variant Calling and Annotation

This repository provides a **step-by-step pipeline** for variant discovery and functional annotation using **GATK** and **ANNOVAR**.  
The workflow starts from reference genome setup and raw reads QC, and ends with **annotated variants** ready for downstream analysis.

---

## 📂 Repository Structure

- `data/` –    (FASTA, BAM, VCF)  
- `scripts/` – Numbered shell scripts (01–14) for each stage of the pipeline  
- `results/` – (Annotated VCF)  
- `README.md` – Project overview and usage instructions  

---

## 🚀 Workflow Overview

### 1. Reference Preparation
- **01_get_ref_seq.sh** → Download human reference genome (hg38).  
- **02_indexing_ref_seq.sh** → Index reference genome with `samtools`.  
- **03_create_dict.sh** → Create sequence dictionary for GATK.  
- **05_bwa_indexing.sh** → Index reference for read alignment.  

### 2. Quality Control
- **04_fastqc.sh** → Run `FastQC` on raw sequencing reads.  

### 3. Alignment & Preprocessing
- **06_mapping_tobamsorting.sh** → Map reads to reference (BWA-MEM) and sort BAM.  
- **07_mark_duplicates.sh** → Mark PCR duplicates.  

### 4. Base Quality Score Recalibration (BQSR)
- **08_get_known_sites_vcf.sh** → Download known sites (dbSNP, Mills, etc.).  
- **09_subsetingforch17.sh** → Subset known sites for chromosome 17 (example).  
- **10_run_bsqr.sh** → Apply BQSR with GATK.  

### 5. Variant Calling
- **11_run_hapotype.sh** → Call variants with `HaplotypeCaller` (gVCF mode).  
- **12_genotyping.sh** → Perform joint genotyping (`GenotypeGVCFs`).  

### 6. Variant Filtering
- **13_hard_filter.sh** → Apply hard filters to SNPs and indels separately, then merge.  

### 7. Annotation
- **14_Annovar_annotation.sh** → Annotate final filtered variants using **ANNOVAR** with `refGene`, ClinVar, dbSNP, etc.  

---

## ⚙️ Requirements

### Software
- [GATK ≥ 4.x](https://gatk.broadinstitute.org/)  
- [ANNOVAR](https://annovar.openbioinformatics.org/)  
- [samtools](http://www.htslib.org/)  
- [BWA](http://bio-bwa.sourceforge.net/)  
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  

### Input Data
- Human reference genome (hg38 FASTA + indices)  
- Known sites VCFs (dbSNP, Mills indels)  
- Raw sequencing reads (FASTQ) or aligned BAM  

### Output
- Raw and filtered VCFs  
- Annotated variant tables (`.hg38_multianno.txt`)  

---

## 📜 License

This project is released for **educational and research purposes**.  
Please check the license terms of **GATK** and **ANNOVAR** before commercial use.  

---

## ✨ Citation

If you use this pipeline, please cite:  
- **McKenna et al., 2010, Genome Analysis Toolkit**  
- **Wang et al., 2010, ANNOVAR: functional annotation of genetic variants**  
