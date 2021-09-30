# Usage:
# Rscript --vanilla generateCache.R -o /mnt/onco0002/repository/software/ensembl-vep/Plugins/

# Argument parser -------------------------------------------------------------------------------------------------

args <- optparse::parse_args(
    optparse::OptionParser(
        add_help_option = T, 
        usage = 'Download annotation files and prepare them for use in VEP.',
        option_list = list(
            optparse::make_option(c('-o', '--outputFolder'), default = NULL, type = 'character', help = 'Output folder', metavar = 'character'),
            optparse::make_option(c('-r', '--removeRawFiles'), default = FALSE, type = 'character', help = 'Remove large raw files (gnoMAD)?', metavar = 'logical')
        )
    )
)

if(is.null(args$outputFolder)) stop('-o / --outputFolder is required.')
if(Sys.which('bcftools') == '') stop('Please install bcftools.')
if(Sys.which('bgzip') == '') stop('Please install bgzip')
if(Sys.which('tabix') == '') stop('Please install tabix')


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
    'gnoMAD - Exome Index' = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi',
    'gnoMAD - Genome' = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz',
    'gnoMAD - Genome Index' = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi',
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
    'gnoMAD - Exome Index' = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi',
    'gnoMAD - Genome' = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz',
    'gnoMAD - Genome Index' = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi',
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
invisible(sapply(filesToDownload$GRCh38, function(x) downloadFile(x, file.path(args$outputFolder, 'GRCh38'))))


# Prepare GENCODE -------------------------------------------------------------------------------------------------

ParallelLogger::logInfo('\t- Preparing GENCODE.')

cleanGENCODE <- function(x){
    
    # Import.
    y <- rtracklayer::import.gff(x)
    
    # Clean-up identifiers (remove ENSEMBL version number) to make them identical to VEP cache.
    y$gene_id <- gsub('\\..*','', y$gene_id)
    y$transcript_id <- gsub('\\..*','', y$transcript_id)
    y$protein_id <- gsub('\\..*','', y$protein_id)
    y$exon_id <- gsub('\\..*','', y$exon_id)
    
    # Remove chr-prefix.
    GenomeInfoDb::seqlevelsStyle(y) <- 'NCBI'
    
    # Export.
    newFile <- gsub(basename(x), paste0('noChrPrefix_', basename(x)), x)
    rtracklayer::export.gff(y, con = newFile, index = T)
    
}

cleanGENCODE(file.path(args$outputFolder, 'GRCh37', basename(filesToDownload$GRCh37$GENCODE)))
cleanGENCODE(file.path(args$outputFolder, 'GRCh38', basename(filesToDownload$GRCh38$GENCODE)))


# # Prepare gnoMAD --------------------------------------------------------------------------------------------------

ParallelLogger::logInfo('\t- Preparing gnoMAD.')

cleangnoMAD <- function(x){
    
    # Perform pre-filtering of PASS-only variant and only retain the AF column.
    cmd <- paste("bcftools view -f PASS", x, "| bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\tAF=%INFO/AF\\n' >> ", gsub('\\.vcf.bgz','_AF.vcf', x))
    base::system(cmd)
    
    # Compress and tabix.
    cmd <- paste('bgzip -c', gsub('\\.vcf.bgz','_AF.vcf', x))
    base::system(cmd)
    cmd <- paste('tabix -p vcf', gsub('\\.vcf.bgz','_AF.vcf.gz', x))
    base::system(cmd)
    
}

cleangnoMAD(file.path(args$outputFolder, 'GRCh37', basename(filesToDownload$GRCh37$`gnoMAD - Exome`)))
cleangnoMAD(file.path(args$outputFolder, 'GRCh38', basename(filesToDownload$GRCh38$`gnoMAD - Exome`)))

cleangnoMAD(file.path(args$outputFolder, 'GRCh37', basename(filesToDownload$GRCh37$`gnoMAD - Genome`)))
cleangnoMAD(file.path(args$outputFolder, 'GRCh38', basename(filesToDownload$GRCh38$`gnoMAD - Genome`)))


# Remove raw files ------------------------------------------------------------------------------------------------

if(args$removeRawFiles == TRUE){
    
    base::system(paste('rm -f', file.path(args$outputFolder, 'GRCh37', basename(filesToDownload$GRCh37$`gnoMAD - Exome`))))
    base::system(paste('rm -f', file.path(args$outputFolder, 'GRCh37', basename(filesToDownload$GRCh37$`gnoMAD - Exome Index`))))
    base::system(paste('rm -f', file.path(args$outputFolder, 'GRCh37', basename(filesToDownload$GRCh37$`gnoMAD - Genome`))))
    base::system(paste('rm -f', file.path(args$outputFolder, 'GRCh37', basename(filesToDownload$GRCh37$`gnoMAD - Genome Index`))))
    
    base::system(paste('rm -f', file.path(args$outputFolder, 'GRCh38', basename(filesToDownload$GRCh37$`gnoMAD - Exome`))))
    base::system(paste('rm -f', file.path(args$outputFolder, 'GRCh38', basename(filesToDownload$GRCh37$`gnoMAD - Exome Index`))))
    base::system(paste('rm -f', file.path(args$outputFolder, 'GRCh38', basename(filesToDownload$GRCh37$`gnoMAD - Genome`))))
    base::system(paste('rm -f', file.path(args$outputFolder, 'GRCh38', basename(filesToDownload$GRCh37$`gnoMAD - Genome Index`))))
    
}

# Finalize --------------------------------------------------------------------------------------------------------

ParallelLogger::logInfo('--- All done! ---')
