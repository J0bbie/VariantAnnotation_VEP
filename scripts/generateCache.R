# Usage:
# Rscript --vanilla generateCache.R -o ~/test

# Argument parser -------------------------------------------------------------------------------------------------

args <- optparse::parse_args(
    optparse::OptionParser(
        add_help_option = T, 
        usage = 'Download annotation files and prepare them for use in VEP.',
        option_list = list(
            optparse::make_option(c('-o', '--outputFolder'), default = NULL, type = 'character', help = 'Output folder', metavar = 'character')
        )
    )
)

if(is.null(args$outputFolder)) stop('-o / --outputFolder is required.')


# Initialize logger -----------------------------------------------------------------------------------------------

ParallelLogger::clearLoggers()

ParallelLogger::registerLogger(
    ParallelLogger::createLogger(
        threshold = "INFO",
        appenders = list(ParallelLogger::createConsoleAppender(layout = ParallelLogger::layoutTimestamp)))
)

ParallelLogger::logInfo('--- Downloading and preparing VEP annotations. ---')

# Generate output folder(s) ---------------------------------------------------------------------------------------

ParallelLogger::logInfo('\t- Generating folders.')

dir.create(file.path(args$outputFolder, 'GRCh37'), recursive = T)
dir.create(file.path(args$outputFolder, 'GRCh38'), recursive = T)


# Required files --------------------------------------------------------------------------------------------------

filesToDownload <- list()

filesToDownload$GRCh37 <- list(
    'GENCODE' = 'http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh37_mapping/gencode.v38lift37.annotation.gtf.gz',
    'CADD - SNV' = 'https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz',
    'CADD - SNV Index' = 'https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz.tbi',
    'CADD - InDel'= 'https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/InDels.tsv.gz',
    'CADD - InDel Index'= 'https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/InDels.tsv.gz.tbi',
    'gnoMAD - Exome' = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz',
    'gnoMAD - Genome' = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz',
    'clinVar' = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20210927.vcf.gz',
    'clinVar Index' = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20210927.vcf.gz.tbi'
)

filesToDownload$GRCh38 <- list(
    'GENCODE' = 'http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v38.annotation.gff3.gz',
    'CADD - SNV' = 'https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz',
    'CADD - SNV Index' = 'https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi',
    'CADD - InDel'= 'https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/InDels.tsv.gz',
    'CADD - InDel Index'= 'https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/InDels.tsv.gz.tbi',
    'gnoMAD - Exome' = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz',
    'gnoMAD - Genome' = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz',
    'clinVar' = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20210927.vcf.gz',
    'clinVar Index' = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20210927.vcf.gz.tbi'
)


# Download files --------------------------------------------------------------------------------------------------

downloadFile <- function(url, outputFolder){
    
    outFile = file.path(outputFolder, basename(url))
    
    if(file.exists(outFile)){
        ParallelLogger::logInfo(sprintf('Already exists, skipping: %s', outFile))
    }else{
        ParallelLogger::logInfo(sprintf('Downloading %s to %s', url, outFile))
        utils::download.file(url, destfile = outFile, method = 'curl', mode = 'wb')
    }
}

ParallelLogger::logInfo('\t- Download files.')

invisible(sapply(filesToDownload$GRCh37, function(x) downloadFile(x, file.path(args$outputFolder, 'GRCh37'))))
invisible(sapply(filesToDownload$GRCh37, function(x) downloadFile(x, file.path(args$outputFolder, 'GRCh38'))))


# Prepare GENCODE -------------------------------------------------------------------------------------------------

ParallelLogger::logInfo('\t- Preparing GENCODE.')

# # Remove chr-prefix to be in concordance with CPCT-02 / HMF VCF files.
# gunzip gencode.v35lift37.annotation.gtf.gz
# grep -v "#" gencode.v35lift37.annotation.gtf | sort -k1,1 -k4,4n -k5,5n -t$'\t' | sed 's/^chr//' | sed 's/^M/MT/' | bgzip -c > noChrPrefix_gencode.v35lift37.annotation.gtf.gz

# 
# z <- rtracklayer::import.gff('/mnt/data/ccbc_environment/general/annotation/hg19/GENCODE/noChrPrefix_gencode.v35lift37.annotation.gtf.gz')
# 
# # Clean-up identifiers (remove ENSEMBL version number) to make them identical to VEP cache.
# z$gene_id <- gsub('\\..*','', z$gene_id)
# z$transcript_id <- gsub('\\..*','', z$transcript_id)
# z$protein_id <- gsub('\\..*','', z$protein_id)
# z$exon_id <- gsub('\\..*','', z$exon_id)
# 
# rtracklayer::export.gff(z, con = '/mnt/data/ccbc_environment/general/annotation/hg19/GENCODE/noChrPrefix_VEP_gencode.v35lift37.annotation.gtf')
# 
# bgzip -c noChrPrefix_VEP_gencode.v35lift37.annotation.gtf > noChrPrefix_VEP_gencode.v35lift37.annotation.gtf.gz
# tabix -p gff noChrPrefix_VEP_gencode.v35lift37.annotation.gtf.gz
# 
# 
# # Prepare gnoMAD --------------------------------------------------------------------------------------------------
 
ParallelLogger::logInfo('\t- Preparing gnoMAD.')

# # # Only keep the max. AF fields and variants passing the gnoMAD filters.
# # bcftools view -h gnomad.exomes.r2.1.1.sites.vcf.bgz > gnomad.exomes.r2.1.1.sites.AF.vcf
# # bcftools view -h gnomad.genomes.r2.1.1.sites.vcf.bgz > gnomad.genomes.r2.1.1.sites.AF.vcf
# # 
# # bcftools view -f PASS gnomad.exomes.r2.1.1.sites.vcf.bgz | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\tAF=%INFO/AF\n' >> gnomad.exomes.r2.1.1.sites.AF.vcf
# # 
# # bcftools view -f PASS gnomad.genomes.r2.1.1.sites.vcf.bgz | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\tAF=%INFO/AF\n' >> gnomad.genomes.r2.1.1.sites.AF.vcf
# # 
# # # Sort and tabix BED file.
# # bgzip gnomad.exomes.r2.1.1.sites.AF.vcf
# # tabix -p vcf gnomad.exomes.r2.1.1.sites.AF.vcf.gz
# # 
# # bgzip gnomad.genomes.r2.1.1.sites.AF.vcf
# # tabix -p vcf gnomad.genomes.r2.1.1.sites.AF.vcf.gz
# # 
# # # Remove original files.
# # rm -f gnomad.exomes.r2.1.1.sites.vcf.bgz* gnomad.genomes.r2.1.1.sites.vcf.bgz*
