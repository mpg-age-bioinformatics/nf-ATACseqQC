#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """

    if [[ "${params.run_type}" == "r2d2" ]] || [[ "${params.run_type}" == "raven" ]] || [[ "${params.run_type}" == "studio" ]]; 
      then
        cd ${params.image_folder}
        if [[ ! -f atacseqqc-1.26.0.sif ]] ;
          then
            singularity pull atacseqqc-1.26.0.sif docker://index.docker.io/mpgagebioinformatics/atacseqqc:1.26.0
        fi
    fi


    if [[ "${params.run_type}" == "local" ]] ; 
      then
        docker pull mpgagebioinformatics/atacseqqc:1.26.0
    fi

    """

}

process atacseqqc_R {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.project_folder}/tmp/ATACSeqQC.${sample}.Rdata").exists() ) 
  input:
    val sample
  
  script:
    """
#!/usr/bin/Rscript

library(openxlsx)
library(ATACseqQC)
library(Rsamtools)
library(ChIPpeakAnno)

if (!dir.exists('/workdir/ATACSeqQC_output')) {
  dir.create('/workdir/ATACSeqQC_output')
}

if (!dir.exists('/workdir/tmp')) {
  dir.create('/workdir/tmp')
}

setwd("/workdir/ATACSeqQC_output/")

# execute ATACSEQQC only if TxDb exists
if("${params.TxDb}" != ""){
library(${params.TxDb})
library(${params.BSgenome})

# library complexity
## input the bamFile from the ATACseqQC package 
bamfile <- '/workdir/bowtie2_output/${sample}.ss.bam'
bamfile.label <- '${sample}'
pdf('/workdir/ATACSeqQC_output/ATACSeqQC.${sample}.pdf')
estimateLibComplexity(readsDupFreq(bamfile))
# fragmentsize distribution
bamfile <- '/workdir/bowtie2_output/${sample}.md.bam'
fragSize <- fragSizeDist(bamfile, bamfile.label)
fragsize = do.call(cbind, fragSize)
write.table(fragsize, '/workdir/ATACSeqQC_output/fragment_size_distribution.${sample}.tsv', sep = '\t', quote = F, row.names = TRUE)
write.xlsx(fragsize, '/workdir/ATACSeqQC_output/fragment_size_distribution.${sample}.xlsx', row.names = TRUE)
## bamfile tags to be read in
possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                 "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                               "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                               "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                               "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                               "U2"))
bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),
                     param = ScanBamParam(tag=unlist(possibleTag)))[[1]][['tag']]
tags <- names(bamTop100)[lengths(bamTop100)>0]
tags
## files will be output into outPath
outPath <- "/workdir/ATACSeqQC_output/splited"
dir.create(outPath)
## shift the coordinates of 5 ends of alignments in the bam file
seqlevelsStyle(${params.atac_genome}) <- "NCBI"
which <- as(seqinfo(${params.atac_genome})[row.names(as.data.frame(seqinfo(${params.atac_genome})))[${params.chrom}]], "GRanges")
shifted_bamfiles = list()
gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
shiftedBamfile <- file.path(outPath, paste0(bamfile.label, ".shifted.bam"))
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)
txdb <- ${params.TxDb}
head(seqlevels(txdb))
seqlevelsStyle(txdb) <- "NCBI"
seqlevels(txdb) <- row.names(as.data.frame(seqinfo(${params.atac_genome})))[${params.chrom}]
txs <- transcripts(txdb)
seqlevels(txs)
# TSS score
tsse <- TSSEscore(gal1, txs)
tsse\$TSSEscore
plot(100*(-9:10-.5), tsse[['values']], type="b",
     xlab="distance to TSS",
     ylab="aggregate TSS score",
     main = '${sample}', pch = 19, cex = 0.7, col = '#87878750')
abline(v = 0, lty = 2, lwd = 0.5)
write.table(as.data.frame(tsse), paste0('/workdir/ATACSeqQC_output/tsse_scores.${sample}.tsv'), row.names = FALSE, sep = '\t', quote = FALSE)
# split bamfile to nucleosome free, mono-nucleosome etc fragments
genome <- ${params.atac_genome}
outPath <- paste0("/workdir/ATACSeqQC_output/splited/${sample}/")
    
#if outputfolder exists, remove and create new
if(file.exists(outPath)){
    # get all files in the directories, recursively
    f <- list.files(outPath, include.dirs = F, full.names = T, recursive = T)
    # remove the files
    print(paste('deleting: ', f))
    file.remove(f)
  }
dir.create(outPath)
  
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = outPath)
# get promoter annotation
setwd(outPath)
library(ChIPpeakAnno)
bamfiles <- c("NucleosomeFree.bam",
              "mononucleosome.bam",
              "dinucleosome.bam",
              "trinucleosome.bam")
TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)
## estimate the library size for normalization
librarySize <- estLibSize(bamfiles)
# not sure what this actually means
cumulativePercentage(bamfiles[1:2], as(seqinfo(${params.atac_genome})[names(which)], "GRanges"))
NTILE <- 101
dws <- ups <- 1010
seqlev = seqlevels(txs)[${params.chrom}]
# calculate density over TSS for all chromosomes
sigs <- enrichedFragments(bamfiles=bamfiles, 
                          TSS=TSS,
                          librarySize=librarySize,
                          seqlev=seqlev,
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)
## log2 transformed signals
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
## get signals normalized for nucleosome-free and nucleosome-bound regions.
out <- featureAlignedDistribution(sigs, 
                                  reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l", 
                                  ylab="Averaged coverage")
                                  
## rescale the nucleosome-free and nucleosome signals to 0~1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)
matplot(out, type="l", xaxt="n", 
        xlab="Position (bp)", 
        ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1, 
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
legend('topleft', c('nucleosome free', 'mononucleosome'), col = c('black', 'red'), lty = c(1,2))
dev.off()
# make TSS heatmap for each chromosome (all chromosomes at once doesn't work; doesn't generate a plot)
# do not do it for the Y chromosome; not enough data
NTILE <- 101
dws <- ups <- 1010
seqlev = seqlevels(txs)[${params.chrom}]
# heatmap folder per sample
dir.create("/workdir/ATACSeqQC_output/chromosome_heatmaps.${sample}/")
for(i in seqlev){
  if( ! i %in% c("Y", "MT")){
    sigs <- enrichedFragments(bamfiles=bamfiles, 
                              TSS=TSS,
                              librarySize=librarySize,
                              seqlev=i,
                              TSS.filter=0.5,
                              n.tile = NTILE,
                              upstream = ups,
                              downstream = dws)
    ## log2 transformed signals
    sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
    #plot heatmap
    pdf(paste0('/workdir/ATACSeqQC_output/chromosome_heatmaps.${sample}/chr_', i, '.pdf'))
    featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                          zeroAt=.5, n.tile=NTILE)
    dev.off()    
    
  }
} 
} else {
  print("THERE IS NO ANNOTATION DATABASE FOR ${params.organism}")
}
save.image("/workdir/tmp/ATACSeqQC.${sample}.Rdata")

tar("chromosome_heatmaps.${sample}.tar.gz", files = "chromosome_heatmaps.${sample}", compression = "gzip")

sessionInfo()

    """
}




workflow images {
  main:
    get_images()
}

workflow {
  if (params.seq == "paired" && params.TxDb != "" && params.BSgenome != "") 
  {
    sample=Channel.fromPath( "${params.project_folder}/bowtie2_output/*.md.bam" )

    sample=sample.map{ "$it.baseName" }
    sample=sample.map{ "$it".replace(".bam","").replace(".md","") }
    sample.view()

    atacseqqc_R(sample)
  }
}

