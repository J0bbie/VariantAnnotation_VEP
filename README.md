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

The following command can be used if all requirements are met:

```bash
# Install cache for GRCh37
perl INSTALL.pl --NO_TEST --NO_HTSLIB --AUTO alcf --PLUGINS SingleLetterAA --CACHEDIR /mnt/onco0002/repository/general/annotation/VEP/GRCh37 --PLUGINSDIR /mnt/onco0002/repository/software/ensembl-vep/Plugins/GRCh37/ --CONVERT --SPECIES homo_sapiens_vep_104_GRCh37

# Install cache for GRCh38, this needs to be a separate folder.
perl INSTALL.pl --NO_TEST --NO_HTSLIB --AUTO alcf --PLUGINS SingleLetterAA --CACHEDIR /mnt/onco0002/repository/general/annotation/VEP/GRCh38 --PLUGINSDIR /mnt/onco0002/repository/software/ensembl-vep/Plugins/GRCh38/ --CONVERT --SPECIES homo_sapiens_vep_104_GRCh37


```

### **Configuration of additional annotations**

Additional files will be downloaded and further processed using `scripts/generateCache.R`. 
This will generate all required files in a user-defined folder (cache) which can be used during annotation.

Currently, the following (non-standard) annotations are added:

- **GENCODE**
  - The default GENCODE version of VEP (GRCh37) is kept at v19. Hence, we utilize a VEP-friendly custom GTF (latest version of GENCODE) to overwrite the default annotations in case of discrepancy.
- **gnoMAD**
  - To add the gnoMAD Allele Frequencies, the gnoMAD VCFs (SNV and InDels) will be processed to greatly reduce their size. Only the genomic positions, REF, ALT and AF will be retained.
- **ClinVar**
  - Required the ClinVar VCF files.
  
## **Generate the VEP command(s)**

Using `scripts/performVEP.R`, we can generate the corresponding VEP command for either GRCh37 or GRCh38 with the respective (additional) annotations.

Example command (GRCh37):
```bash
Rscript --vanilla scripts/performVEP.R -b GRCh37 -i  ~/test/CPCT02010257T.purple.somatic.vcf.gz -x /mnt/onco0002/repository/software/ensembl-vep/vep -g /mnt/onco0002/repository/software/ensembl-vep/Plugins/GRCh37/noChrPrefix_gencode.v38lift37.annotation.gtf.bgz -p /mnt/onco0002/repository/software/ensembl-vep/Plugins/ -c /mnt/onco0002/repository/general/annotation/VEP/ -f /mnt/onco0002/repository/general/annotation/VEP/GRCh37/homo_sapiens/104_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
```

Example command (GRCh38):
```bash
Rscript --vanilla scripts/performVEP.R -b GRCh38 -i asd.vcf -x /mnt/onco0002/repository/software/ensembl-vep/vep -g /mnt/onco0002/repository/software/ensembl-vep/Plugins/GRCh38/noChrPrefix_gencode.v38.annotation.gff3.bgz -p /mnt/onco0002/repository/software/ensembl-vep/Plugins/ -c /mnt/onco0002/repository/general/annotation/VEP/ -f /mnt/onco0002/repository/general/annotation/VEP/GRCh38/homo_sapiens/104_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
```
