#### poisson
args=(commandArgs(TRUE))


source("./library_functions_ZCL_chip.R")


Bed<-args[1]      #Bed file normalization
ALL_PEAKS<-args[2]   #TRUE or FALSE
klustering<-args[3]     #numeric, from 1 to whatever
control<-args[4]  #TRUE or FALSE
anchura_max<-args[5]
anchura_min<-args[6]
cat(Bed,"\n")
cat(ALL_PEAKS,"\n")
 cat(klustering, "\n")
 cat(control,"\n")
cat(anchura_min, "\n")
cat(anchura_max,"\n")
peaks_clustered<-read.table(Bed, quote="", sep="\t")
#peaks_clustered$width<-peaks_clustered$V3-peaks_clustered$V2
#peaks_clustered<-peaks_clustered[which(peaks_clustered$width < as.numeric(paste(anchura_max))),]
#peaks_clustered<-peaks_clustered[which(peaks_clustered$width > as.numeric(paste(anchura_min))),]
peaks_clustered$poisson<-1-(ppois(as.numeric(paste(peaks_clustered[,5])), as.numeric(paste(peaks_clustered[,6]))))
peaks_clustered$poisson_FDR<-p.adjust(peaks_clustered$poisson)

if(ALL_PEAKS==FALSE){
    selected<-peaks_clustered[which(peaks_clustered$poisson_FDR < 0.01),]
    selected$width<-selected[,3] - selected[,2]
    #selected<- selected[which(selected$width < as.numeric(paste(anchura_max))),]
    output<-data.frame("chr"=NULL,"start"=NULL,"end"=NULL, "area"=NULL, "area_input"=NULL)
    for (chr in unique(paste(selected$V1))){
        cat(chr,"\r")
        flush.console()
        chunk_df<-selected[which(selected$V1==paste(chr)),]
        if(nrow(chunk_df)>0){
            chunk_df_ranges<-IRanges(start=chunk_df$V2,end=chunk_df$V3)
            chunk_df_clustered<-reduce(chunk_df_ranges, min.gapwidth=as.numeric(paste(klustering)))
            chunk_df_clustered<-as.data.frame(chunk_df_clustered)
            chunk_df_clustered<-unique(chunk_df_clustered)
            chr_df<-data.frame("V1"= rep(chr, nrow(chunk_df_clustered)))
            output_chr<-cbind(chr_df, chunk_df_clustered)
            file_dir<-dirname(Bed)
            file_dir<-dirname(file_dir)
            file_dir<-paste0(file_dir,"/TMP/chip/")
            file_dir_less<-dirname(file_dir)
            file_dir_input<-paste0(file_dir_less,"/input/")
            file_chip<-paste0(file_dir, chr,".chrom")
            chipchrom<-fread(paste(file_chip))
            chipchrom<-as.data.frame(chipchrom)
            quantified<-quantification(output_chr,chipchrom$V3)
            output_chr$area_chip<-log2(unlist(quantified))

            if(control==TRUE){
            file_input<-paste0(file_dir_input, chr,".chrom")
            factor_file<-paste0(file_dir, chr, ".fc")
            factor<-read.table(factor_file)
            fc<-as.numeric(paste(factor$V1))
            controlchrom<-fread(paste(file_input))
            controlchrom<-as.data.frame(controlchrom)
            quantified_input<-quantification(output_chr,controlchrom$V3)
            output_chr$area_input<-log2(unlist(quantified_input)) +log10(fc)
            output_chr$area_input[!is.finite(output_chr$area_input)]<-0

            output<-rbind(output, output_chr)
            }
        }
    }

    output$name<-paste0("peak_",seq(1:nrow(output)),"@",output$area_chip,"_",output$area_input)
    output_prepared<-output[,c("V1","start","end","name")]
    output_prepared$width<-output_prepared$end - output_prepared$start
    output_prepared<-output_prepared[which(output_prepared$width < as.numeric(paste(anchura_max))),]
    output_prepared<-output_prepared[which(output_prepared$width > as.numeric(paste(anchura_min))),]

    write.table(output_prepared, file=paste0(lapply(strsplit(paste(Bed),".model"),"[",1),"_output_poissonFDR.bed"),col.names=F, row.names=F, quote=F, sep="\t")

}else{

    output<-data.frame("chr"=NULL,"start"=NULL,"end"=NULL, "area"=NULL)
    for (chr in unique(paste(peaks_ranges$V1))){
        chunk_df<-peaks_clustered[which(peaks_clustered$V1==paste(chr)),]
        chunk_df_ranges<-IRanges(start=chunk_df$V2,end=chunk_df$V3)
        chunk_df_clustered<-reduce(chunk_df_ranges, min.gapwidth=as.numeric(paste(klustering)))
        chunk_df_clustered<-as.data.frame(chunk_df_clustered)
        chunk_df_clustered<-unique(chunk_df_clustered)
        chr_df<-data.frame("V1"= rep(chr, nrow(chunk_df_clustered)))
        output_chr<-cbind(chr_df, chunk_df_clustered)
        file_dir<-dirname(Bed)
        file_dir<-dirname(file_dir)
        file_dir<-paste0(file_dir,"/TMP/chip/")
        file_chip<-paste0(file_dir, chr,".chrom")
        chipchrom<-fread(paste(file_chip))
        chipchrom<-as.data.frame(chipchrom)
        quantified<-quantification(output_chr,chipchrom$V3)
        output_chr$area<-log2(unlist(quantified))
        output<-rbind(output, output_chr)
    }
    output$name<-paste0("peak_",seq(1:nrow(output)),"@",output$area_chip)
    output_prepared<-output[,c("chr","start","end","name")]
    output_prepared$width<-output_prepared$end - output_prepared$start
    #output_prepared<-output_prepared[which(output_prepared$width < as.numeric(paste(anchura_max))),]
    #output_prepared<-output_prepared[which(output_prepared$width > as.numeric(paste(anchura_min))),]

    write.table(output_prepared, file=paste0(lapply(strsplit(paste(Bed),".model"),"[",1),"_output_poissonFDR.bed"),col.names=F, row.names=F, quote=F, sep="\t")

}
