# Workflow - Annotation using VEP
![license](https://img.shields.io/badge/license-GPL--3-blue.svg)

![Maintainer](https://img.shields.io/badge/Job%20van%20Riet-lightgrey.svg?label=Maintainer)


### **Table of Content**

- [Workflow - Annotation using VEP](#workflow---annotation-using-vep)
    - [**Table of Content**](#table-of-content)
  - [**Introduction**](#introduction)
  - [**Installation of VEP**](#installation-of-vep)
  - [**Configuration of additional annotations**](#configuration-of-additional-annotations)


## **Introduction**
Workflow for annotating (somatic) variants using [Variant Effect Predictor (VEP)](https://github.com/Ensembl/ensembl-vep) for GRCh37 or GRCh38.

Briefly, it will perform the following annotations:
- Standard annotations (using custom GTF for GRCh37)
- Combined Annotation Dependent Depletion (CADD)
- gnoMAD Allele Frequencies - Genome and Exome (v2.1.1)
- ClinVar
- SingleLetterAA

The scripts within this workflow will generate a folder (cache) containing all the annotation files required for VEP and additional downstream analysis.

## **Installation of VEP**

Following the [author's instructions](https://github.com/Ensembl/ensembl-vep), install VEP and required CPAN plugins.

In addition, the following tools are also required:
- bcftools
- bgzip
- tabix

## **Configuration of additional annotations**

Additional files will be downloaded and further processed using `scripts/generateCache.R`. This will generate all required files in a user-defined folder (cache) which can be used in annotation.

- **GENCODE**
  - The default GENCODE version of VEP (GRCh37) is kept at v19. Hence, we utilize a VEP-friendly custom GTF (latest version of GENCODE) to overwrite the default annotations in case of discrepancy.
- **CADD**
  - The Combined Annotation Dependent Depletion (CADD) requires two additional files (per reference) containing the raw and phred scores for each SNV and InDel.
- **gnoMAD**
  - To add the gnoMAD Allele Frequencies, the gnoMAD VCFs (SNV and InDels) will be processed to greatly reduce their size. Only the genomic positions, REF, ALT and AF will be retained.
- **ClinVar**
  - Required the ClinVar VCF files.
