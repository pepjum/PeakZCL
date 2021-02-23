#!/usr/bin/env Rscript

#source("/home/bioinformatica/datos/03_Analysis/jgonzalez.69/scripts_chipseq/library_functions_ZCL_chip.R")
source("./library_functions_ZCL_chip.R")

#change this path for your own

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="PATH to a chipseq file name", metavar="character"),
  make_option(c("-c", "--control"), type="character", default=NULL,
              help="PATH to a control file name", metavar="character"),
  make_option(c("-n", "--normalization"), type="character", default="SES",
              help="Select the signal normalization method.Options are SES,SDS an N (without).Default method is SES", metavar="character"),
  make_option(c("-w", "--Decimating"), type="numeric", default=50,
              help="Decimating signal factor. If you decide not to decimate the original signal, the program will take a long time. If the value is high, you can lose a lot of information. Default is 50", metavar="character"),
  make_option(c("-l", "--length_signal_wavelet"), type="numeric", default=200000,
              help="Length signal wavelet.Default value is 2000000.If signal is longer than 4000000, wavelet function will not work.If the value is less than 100000, the peak calling function may cause errors", metavar="character"),
  make_option(c("-s", "--scale"), type="numeric", default=100,
              help="Number of scales of wavelet. Default is 100", metavar="character"),
  make_option(c("-g", "--Levth"), type="numeric", default=5,
              help="Zero-crossing lines leverage threshold. Default is 5", metavar="character"),
  make_option(c("-z", "--enrichment"), type="character", default="",
              help="histone or TF. Default=TF", metavar="character"),
  make_option(c("-j", "--clustering"), type="numeric", default=300,
              help="Bp range between peaks to clustering purposes. Default is 300", metavar="character"),
  make_option(c("-a", "--area_chip"), type="numeric", default=0.69,
              help=" Minimum value of the area of the selected peaks. The value is given in log2(x). Default is 0.69", metavar="character"),
  make_option(c("-k", "--minpeak"), type="numeric", default=10,
              help="minimum width of the peak. Default is 10", metavar="character"),
  make_option(c("-x", "--maxpeak"), type="numeric", default=500000,
              help="minimum width of the peak. Default is 500000", metavar="character"),
  make_option(c("-o", "--out_file"), type="character",default=NULL,
              help="Output BED file name.Default name will be the same at your experiment.", metavar="character"),
  make_option(c("-m", "--all_peaks"), type="character", default=FALSE,
              help="Show all peaks detected", metavar="character")

);

threshold<-4

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# cat("file", opt$file, "\n")
# cat("control", opt$control, "\n")


if(is.null(opt$file)){
  print_help(opt_parser)
  stop("At least two argument must be supplied (input file and output directory)","\n",
  "Example of use: Rscript (PATH to ZCL_Rscript file in your system) -f (PATH to chipseq file) -c (PATH to control file), -z (PATH to chromsizes) -d (PATH to output directory)", call.=FALSE)
}

fileName_chip<-basename(opt$file)
# cat("loading chip file ",fileName_chip,"\n")
chipchrom<-fread(opt$file)
chipchrom<-as.data.frame(chipchrom)
chr_number<-as.character(unlist(chipchrom$V1[1]))

out_dir<-dirname(opt$file)
#cat("outdir ", out_dir, "\n")
# cat("decimating chip signal...",fileName_chip,"\n")
if(opt$Decimating==0){
    decimated_signal_chip<-chipchrom$V3
}else{
    decimated_signal_chip<-decimate(chipchrom$V3,as.numeric(paste(opt$Decimating)))
    decimated_signal_chip[decimated_signal_chip<0]<-0
}

if(!is.null(opt$control)){
    fileName_control<-basename(opt$control)
    # cat("loading control file ",fileName_control,"\n")
    controlchrom<-fread(opt$control)
    controlchrom<-as.data.frame(controlchrom)

    noiselevel<-1.5*(sd(controlchrom$V3))
    if(opt$enrichment=="histone"){
        if( noiselevel < threshold){
            noiselevel<- noiselevel
        }

    }else{  noiselevel<-noiselevel}
}else{
    cat("estimating noise level from signal","\n")
    # max_index<- nrow(chipchrom)-1000001
    # index_creator<-sample(1:max_index, 50)
    signal_splitted<-split(chipchrom$V3, ceiling(seq_along(chipchrom$V3)/(length(chipchrom$V3)/50)))
    selected_signal<-sample(signal_splitted, 20)

    medias<-c()
    medianas<-c()
    standard_desviation<-c()
    for(i in 1:length(selected_signal)){
        region_signal<-selected_signal[[i]]
        media_signal<-mean(region_signal)
        mediana_signal<-median(region_signal)
        desviation_signal<-sd(region_signal)
        medias<-c(medias, media_signal )
        medianas<-c(medianas, mediana_signal)
        standard_desviation<-c(standard_desviation, desviation_signal)
    }
    delta_m<-c()
    for(i in 1:length(medias)){
        delta<-(medias[i])-(medianas[i])
        delta_m<-c(delta_m, delta)
    }
    cutoff_value<-summary(delta_m)[4]
    desviation_selected<-c()
    for(i in 1:length(medias)){
        if((medias[i])-(medianas[i]) < cutoff_value){
            desviation_selected<-c(desviation_selected,standard_desviation[i])
        }
    }

#    noiselevel<-sum(desviation_selected)/length(desviation_selected)
    noiselevel<-1.5*(sum(desviation_selected)/length(desviation_selected))
    if(opt$enrichment=="histone"){
        if( noiselevel < threshold){
            noiselevel<- noiselevel
        }

    }else{  noiselevel<-noiselevel}

}

# cat("estimating CWT and peak calling...",fileName_chip,"\n")
start_end_peaks<-estimating_CWT_and_peak_calling(decimated_signal_chip,opt$length_signal_wavelet,opt$scale,opt$Levth,opt$Decimating, opt$Threshold, noiselevel)
peaks_in_position<-((start_end_peaks)-1)*as.numeric(paste(opt$Decimating))


#start_end_peaks_ranges<-IRanges(start=start_end_peaks$start,end=start_end_peaks$end)
peaks_ranges<-IRanges(start=peaks_in_position$start,end=peaks_in_position$end)

peaks_clustered<-reduce(peaks_ranges, min.gapwidth=opt$clustering)
peaks_clustered<-as.data.frame(peaks_clustered)
peaks_clustered<-unique(peaks_clustered)


if(!is.null(opt$control)){

  quantification_bins_chip<-list()
  n<-1000
  k<-1
  max_len<-length(chipchrom$V3)
  for(i in seq(from=1, to=length(chipchrom$V3), by=n)){
        if((i)+n > (max_len)){
            bin<-chipchrom$V3[(i):max_len]
            sum_bin<-sum(as.numeric(paste(bin)))
            quantification_bins_chip[[k]]<-sum_bin
            # cat(i,"_",max_len, " / ", max_len,"\r")
            # flush.console()
    }
    else{
            bin<-chipchrom$V3[(i):(i + n-1)]
            sum_bin<-sum(as.numeric(paste(bin)))
            quantification_bins_chip[[k]]<-sum_bin
            k<-k+1
            # cat((i),"_",(i+n-1), " / ", (max_len),"\r")
            # flush.console()
    }

  }

  quantification_bins_control<-list()
  n<-1000
  k<-1
  max_len<-length(controlchrom$V3)
  for(i in seq(from=1, to=length(controlchrom$V3), by=n)){
        if((i)+n > (max_len)){
            bin<-controlchrom$V3[(i):max_len]
            sum_bin<-sum(as.numeric(paste(bin)))
            quantification_bins_control[[k]]<-sum_bin
            # cat(i,"_",max_len, " / ", max_len,"\r")
            # flush.console()
    }
    else{
            bin<-controlchrom$V3[(i):(i +n-1)]
            sum_bin<-sum(as.numeric(paste(bin)))
            quantification_bins_control[[k]]<-sum_bin
            k<-k+1
            # cat((i),"_",(i+n-1), " / ", (max_len),"\r")
            # flush.console()
    }

  }
  quantification_bins_chip<-unlist(quantification_bins_chip)
  quantification_bins_control<-unlist(quantification_bins_control)

  if(opt$normalization=="SES"){
    # cat("normalizing signal....",fileName_control,"with ", opt$normalization,"method","\n")
    factor_x<-NormalizationChIP(as.numeric(paste(quantification_bins_chip)),as.numeric(paste(quantification_bins_control)),method="SES")
  }else if(opt$normalization=="SDS"){
    # cat("normalizing signal....",fileName_control,"with ", opt$normalization,"method","\n")
    factor_x<-NormalizationChIP(as.numeric(paste(quantification_bins_chip)),as.numeric(paste(quantification_bins_control)),method="SDS")
  }else if(opt$normalization=="N"){
    factor_x<-1
  }
  c_normalized<-as.numeric(paste(controlchrom$V3)) *factor_x
  controlchrom$V3<-c_normalized
# normalizacion nueva

   if (opt$all_peaks==FALSE){
    peaks_clustered<-peaks_clustered[which(peaks_clustered$width > opt$minpeak),]
    peaks_clustered<-peaks_clustered[which(peaks_clustered$width < opt$maxpeak),]
    quantified<-quantification(peaks_clustered,chipchrom$V3)
    quantified_input<-quantification(peaks_clustered, controlchrom$V3)
    peaks_clustered$area_chip<-as.numeric(paste(quantified))
    peaks_clustered$area_input<-as.numeric(paste(quantified_input))
    #peaks_clustered<-peaks_clustered[which(peaks_clustered$area_chip > opt$area_chip),]
    }
    bed<-data.frame("chr"=rep(chr_number, nrow(peaks_clustered)), "start"=peaks_clustered$start, "end"= peaks_clustered$end, "name"= paste0("peak_",seq(1:nrow(peaks_clustered)), "@",peaks_clustered$area_chip,"_", peaks_clustered$area_input))

    write.table(bed,file=paste0(out_dir,"/",chr_number,".bed"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(factor_x,file=paste0(out_dir,"/",chr_number,".fc"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}else{
    quantified<-quantification(peaks_clustered,chipchrom$V3)
    peaks_clustered$quantification<-log10(unlist(quantification))
    bed<-data.frame("chr"=rep(chr_number, nrow(peaks_clustered)), "start"=peaks_clustered$start, "end"= peaks_clustered$end, "name"= paste0("peak_",seq(1:nrow(peaks_clustered)),"@",peaks_clustered$quantification))
    write.table(bed,file=paste0(out_dir,"/",chr_number,".bed"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

}

cat("Peak detection of", fileName_chip,"DONE!","\n")
