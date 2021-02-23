#### poisson
source("./library_functions_ZCL_chip.R")


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="peaks chip list in bed format", metavar="character"),
  make_option(c("-c", "--counts_chip"), type="character", default=NULL,
              help="Counts chip", metavar="character"),
  make_option(c("-d", "--counts_input"), type="character", default=NULL,
              help="Counts input file", metavar="character"),
  make_option(c("-n", "--factor"), type="character", default=NULL,
              help="normalization file factor", metavar="character"),
  make_option(c("-p", "--outdir"), type="character", default=NULL,
              help="output dir", metavar="character"),
  make_option(c("-o", "--outnamefile"), type="character", default=NULL,
              help="Output name file", metavar="character"),
  make_option(c("-m", "--allpeaks"), type="character", default=FALSE,
              help="returns all peaks detected or not", metavar="character"),
  make_option(c("-b", "--klustering"), type="numeric", default=1,
              help="range between peaks for clustering purposes", metavar="character"),
  make_option(c("-w", "--chip"), type="character", default=NULL,
              help="chip file coverage", metavar="character")


);



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# cat(opt$outdir,"\n")
# cat(opt$outnamefile,"\n")
# cat(opt$file,"\n")
# cat(opt$counts_chip,"\n")
# cat(opt$counts_input,"\n")
#

if(is.null( opt$counts_input)){
    bedfile<-read.table(opt$file)
    counts_chip<-read.table(opt$counts_chip)
    counts_chip<-counts_chip[,5]
    df_tot<-cbind(bedfile, counts_chip)
    write.table(df_tot, file=paste0(opt$outdir, opt$outnamefile,".Norm"),col.names=F ,row.names=F, quote=F, sep="\t")

}else{
    bedfile<-read.table(opt$file)
    counts_chip<-read.table(opt$counts_chip)
    counts_chip<-counts_chip[,5]
    counts_input<-read.table(opt$counts_input)
    counts_input<-counts_input[,5]
    factor<-read.table(opt$factor)
    factor<-as.numeric(paste(factor$V1))
    counts_input<-(counts_input)*factor

    peaks_clustered<-cbind(bedfile, counts_chip, counts_input)

    peaks_clustered$areas<-lapply(strsplit(paste(peaks_clustered$V4),"@"),"[",2)
    peaks_clustered$area_chip<-as.numeric(paste(lapply(strsplit(paste(peaks_clustered$areas),"_"),"[",1)))
    peaks_clustered$area_input<-as.numeric(paste(lapply(strsplit(paste(peaks_clustered$areas),"_"),"[",2)))
    peaks_clustered$area_chip<-log2(peaks_clustered$area_chip)
    peaks_clustered$area_input<-(peaks_clustered$area_input)*factor

    peaks_clustered$area_input<-log2(peaks_clustered$area_input)
    peaks_clustered$area_input[!is.finite(peaks_clustered$area_input)]<-0
    peaks_clustered$fold_change_log2<-peaks_clustered$area_chip - peaks_clustered$area_input
    peaks_clustered<-peaks_clustered[which(peaks_clustered$fold_change > 0),]
    #peaks_clustered$fold_change_log2<-log2(peaks_clustered$fold_change)
#     peaks_clustered$fold_change_counts<-peaks_clustered$counts_chip - peaks_clustered$counts_input
# #    peaks_clustered<-peaks_clustered[which(peaks_clustered$fold_change > 0),]
#     peaks_clustered$zscore<-scale(peaks_clustered$fold_change)
#     peaks_clustered$zscore_counts<-scale(peaks_clustered$fold_change_counts)
#     #peaks_clustered$zscore_arr<-(peaks_clustered$zscore)+1
#     peaks_clustered$areas<-NULL
#     selected<-peaks_clustered[which(peaks_clustered$zscore > abs(qnorm(0.1))),]
    # selected<-peaks_clustered
    # selected$area_chip<-NULL
    # selected$area_input<-NULL
    # selected$fold_change<-NULL
    # selected$zscore<-NULL
    # #selected$zscore_arr<-NULL

    # filtering by fold_change

    selected<-peaks_clustered
    selected$areas<-NULL

    if(nrow(selected) >0){
    write.table(selected, file=paste0(paste(opt$outdir), paste(opt$outnamefile),".Norm"),col.names=F ,row.names=F, quote=F, sep="\t")
    }else{
        break
    }
}

#### fold change, zscore
