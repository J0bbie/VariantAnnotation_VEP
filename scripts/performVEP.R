# Argument parser -------------------------------------------------------------------------------------------------

args <- optparse::parse_args(
    optparse::OptionParser(
        add_help_option = T,
        usage = 'Perform ENSEMBL VEP using standard cache and plugins on VCFs within a given folder.',
        option_list = list(
            optparse::make_option(c('-b', '--build'), default = 'GRCh37', type = 'character', help = 'Output folder', metavar = 'character'),
            optparse::make_option(c('-o', '--outputFolder'), default = NULL, type = 'character', help = 'Output folder (if left empty, will use the same folder as input VCF)', metavar = 'character'),
            optparse::make_option(c('-t', '--threads'), default = 4, type = 'character', help = 'Nr. of threads per file', metavar = 'character'),
            
            optparse::make_option(c('-x', '--pathVEP'), default = NULL, type = 'character', help = 'Path to VEP', metavar = 'character'),
            optparse::make_option(c('-i', '--pathVCF'), default = NULL, type = 'character', help = 'Path to VCF file', metavar = 'character'),
            optparse::make_option(c('-f', '--pathFASTA'), default = NULL, type = 'character', help = 'Path to genome FASTA file', metavar = 'character'),
            optparse::make_option(c('-g', '--pathGTF'), default = NULL, type = 'character', help = 'Path to GTF file', metavar = 'character'),
            optparse::make_option(c('-c', '--pathCACHE'), default = NULL, type = 'character', help = 'Path to VEP cache folder, e.g. /homo_sapiens/104_GRCh37/', metavar = 'character'),
            optparse::make_option(c('-p', '--pathPLUGIN'), default = NULL, type = 'character', help = 'Path to VEP plugins (and respective data)', metavar = 'character')
        )
    )
)

if(is.null(args$pathVEP)) stop('-x / --pathVEP is required.')
if(is.null(args$pathVCF)) stop('-i / --pathVCF is required.')
if(is.null(args$pathCACHE)) stop('-c / --pathCACHE is required.')
if(is.null(args$pathFASTA)) stop('-f / --pathFASTA is required.')
if(is.null(args$pathGTF)) stop('-g / --pathGTF is required.')
if(is.null(args$pathPLUGIN)) stop('-p / --pathPLUGIN is required.')

# Combine plugin and build.
args$pathPLUGIN <- file.path(args$pathPLUGIN, args$build)
args$pathCACHE <- file.path(args$pathCACHE, args$build)

# Find synonyms file.
args$pathSYNONYM <- list.files(list.files(file.path(args$pathCACHE, 'homo_sapiens'), full.names = T), pattern = 'chr_synonyms.txt', full.names = T, all.files = F, include.dirs = F, recursive = F)

# Determine output file -------------------------------------------------------------------------------------------

outFile <- base::gsub('\\.*vcf.*', '_annotedWithVEP.vcf.gz', base::basename(args$pathVCF))

if(is.null(args$outputFolder)) outFile <- file.path(base::dirname(args$pathVCF), outFile)
if(!is.null(args$outputFolder)) outFile <- file.path(args$outputFolder, outFile)


# Generate VEP command --------------------------------------------------------------------------------------------

cmd.VEP <- base::sprintf('%s -i %s --fasta %s --fork %s --gtf %s --dir_plugins %s --assembly %s --dir_cache %s -o %s --synonyms %s --dir_plugins %s
--custom %s/gnomad.exomes.r2.1.1.sites_AF.vcf.gz,gnomADe,vcf,exact,0,AF 
--custom %s/gnomad.genomes.r2.1.1.sites_AF.vcf.gz,gnomADg,vcf,exact,0,AF 
--custom %s/clinvar_20210927.vcf.gz,ClinVar,vcf,exact,0,CLNDN,CLNHGVS,CLNSIG 
--plugin SingleLetterAA 
--ccds --hgvs --symbol --force_overwrite --numbers --domains --canonical --protein --biotype --uniprot --tsl --appris --gene_phenotype --pubmed --variant_class 
--vcf_info_field ANN --cache --offline --no_stats 
--quiet --format vcf --vcf --compress_output bgzip', 
              args$pathVEP, args$pathVCF, args$pathFASTA, args$threads, args$pathGTF, args$pathPLUGIN, args$build, args$pathCACHE, outFile, args$pathSYNONYM, args$pathPLUGIN, args$pathPLUGIN, args$pathPLUGIN, args$pathPLUGIN)

base::system(paste('echo', gsub('\n', '', cmd.VEP)))