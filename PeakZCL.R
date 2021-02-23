#!/usr/bin/env Rscript

#source("/home/bioinformatica/datos/03_Analysis/jgonzalez.69/scripts_chipseq/library_functions_ZCL_chip.R")
source("./library_functions_ZCL_chip.R")

#change this path for your own

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="PATH to a chipseq file name", metavar="character"),
  make_option(c("-c", "--control"), type="character", default=NULL,
              help="PATH to a control file name", metavar="character"),
  make_option(c("-d", "--output_dir"), type="character", default=NULL,
              help="A directory for ZCL to create its output", metavar="character"),
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
              help=" Minimum value of the area of the selected peaks. The value is given in log10(x). Default is 0.69", metavar="character"),
  make_option(c("-k", "--fold_change"), type="numeric", default=0.1,
              help="Difference between peak area on the chip and control. The value is given in log10(x). Default is 0.1", metavar="character"),
  make_option(c("-o", "--out_file"), type="character",default=NULL,
              help="Output BED file name.Default name will be the same at your experiment.", metavar="character"),
  make_option(c("-m", "--all_peaks"), type="character", default=FALSE,
              help="Show all peaks detected", metavar="character")

);

threshold<-4

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if(is.null(opt$file)){
  print_help(opt_parser)
  stop("At least two argument must be supplied (input file and output directory)","\n",
  "Example of use: Rscript (PATH to ZCL_Rscript file in your system) -f (PATH to chipseq file) -c (PATH to control file), -z (PATH to chromsizes) -d (PATH to output directory)", call.=FALSE)
}

fileName_chip<-basename(opt$file)
cat("loading chip file ",fileName_chip,"\n")
chipchrom<-fread(opt$file)
chipchrom<-as.data.frame(chipchrom)
chr_number<-as.character(unlist(chipchrom$V1[1]))

# cat("decimating chip signal...",fileName_chip,"\n")
if(opt$Decimating==0){
    decimated_signal_chip<-chipchrom$V3
}else{
    decimated_signal_chip<-decimate(chipchrom$V3,as.numeric(paste(opt$Decimating)))
    decimated_signal_chip[decimated_signal_chip<0]<-0
}

if(!is.null(opt$control)){
    fileName_control<-basename(opt$control)
    cat("loading control file ",fileName_control,"\n")
    controlchrom<-fread(opt$control)
    controlchrom<-as.data.frame(controlchrom)
    noiselevel<-sd(controlchrom$V3)
    if(opt$enrichment=="histone"){
        if( noiselevel < threshold){
            noiselevel<- threshold
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
    noiselevel<-sum(desviation_selected)/length(desviation_selected)
    if(opt$enrichment=="histone"){
        if( noiselevel < threshold){
            noiselevel<- threshold
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

signal_selected_chip<-c()
intensities<-c()
#chipchrom$V4<-0
for (i in 1:nrow(peaks_clustered)){
    cat(i, "/", nrow(peaks_clustered),"\r")
    flush.console()
  peaks_signal_ranges_chip<-chipchrom$V3[peaks_clustered$start[i]:peaks_clustered$end[i]]
#  chipchrom$V4[peaks_clustered$start[i]:peaks_clustered$end[i]]<-1
  max_intensities<-max(peaks_signal_ranges_chip)
  intensities<-c(intensities, max_intensities)
  signal_selected_chip<-c(signal_selected_chip,peaks_signal_ranges_chip)
}

if(!is.null(opt$control)){

  signal_selected_input<-c()
 # controlchrom$V4<-0
  for (i in 1:nrow(peaks_clustered)){
      cat(i, "/", nrow(peaks_clustered),"\r")
      flush.console()
      peaks_signal_ranges_input<-controlchrom$V3[peaks_clustered$start[i]:peaks_clustered$end[i]]
#      controlchrom$V4[peaks_clustered$start[i]:peaks_clustered$end[i]]<-1

      signal_selected_input<-c(signal_selected_input,peaks_signal_ranges_input)
      #signal_selected_input[is.na(signal_selected_input)]<-0
  }

  if(opt$normalization=="SES"){
    # cat("normalizing signal....",fileName_control,"with ", opt$normalization,"method","\n")
    factor_x<-NormalizationChIP(as.numeric(paste(chipchrom$V3)),as.numeric(paste(controlchrom$V3)),method="SES")
  }else if(opt$normalization=="SDS"){
    # cat("normalizing signal....",fileName_control,"with ", opt$normalization,"method","\n")
    factor_x<-NormalizationChIP(as.numeric(paste(chipchrom$V3)),as.numeric(paste(controlchrom$V3)),method="SDS")
  }else if(opt$normalization=="N"){
    factor_x<-1
  }

# normalizacion nueva




  c_normalized<-(controlchrom$V3)*as.numeric(paste(factor_x))

  controlchrom$V3<- c_normalized

  intensities_input<-c()
  for (i in 1:nrow(peaks_clustered)){
      peaks_signal_ranges_input<-controlchrom$V3[peaks_clustered$start[i]:peaks_clustered$end[i]]
      intensities_in<-max(peaks_signal_ranges_input)
      intensities_input<-c(intensities_input, intensities_in)

  }

}

# cat("quantification...",fileName_chip,"\n")
quantified<-quantification(peaks_clustered,chipchrom$V3)

if(!is.null(opt$control)){
  quantified_input<-quantification(peaks_clustered,controlchrom$V3)
}

# write.table(cbind(data.frame("chr"=rep("chr3", nrow(peaks_clustered))),peaks_clustered[,c("start","end")], data.frame("name"=seq(1:nrow(peaks_clustered)))), file="chr3_regions.bed", col.names=F, row.names=F, quote=F, sep="\t")
#
# counts_chip<-read.table("~/data/pepe/89_EGR1_pruebas_chr3/chr3_counts_CHIP_multicov.txt")
# counts_input<-read.table("~/data/pepe/89_EGR1_pruebas_chr3/chr3_counts_INPUT_multicov.txt")
#
# peaks_clustered$counts_chip<-as.numeric(paste(counts_chip$V5))
# peaks_clustered$counts_input<-as.numeric(paste(counts_input$V5))
# peaks_clustered$counts_input_N<-peaks_clustered$counts_input)*factor_x
# peaks_clustered$poisson<-1-(ppois(as.numeric(paste(peaks_clustered$counts_chip)), as.numeric(paste(peaks_clustered$counts_input_N))))
# peaks_clustered$poisson_FDR<-p.adjust(peaks_clustered$poisson, method="fdr")
# selected<-peaks_clustered[which(peaks_clustered$poisson_FDR < 0.05),])

# peaks_clustered$diff<-(peaks_clustered$counts_chip) - (peaks_clustered$counts_input)
# peaks_clustered$dpois<-dpois(peaks_clustered$diff, lambda=mean(peaks_clustered$diff))

if(!is.null(opt$control)){

      peaks_clustered$intensity_chip<-unlist(intensities)
      peaks_clustered$intensity_input<-unlist(intensities_input)


      peaks_in_position<-statistics(peaks_clustered,as.numeric(paste(opt$area_chip)),as.numeric(paste(opt$fold_change)),quantified,quantified_input, paste(opt$enrichment))
  # }
}else{
  peaks_in_position<-peaks_clustered
  peaks_in_position$area_chip<-log10(unlist(quantified))
}


options(scipen=100)
bed<-BED_files(peaks_in_position,chr_number,opt$control,opt$all_peaks)

if(nrow(bed) ==0){
    cat("no peaks in ",fileName_chip,"\n")
    break
}

cat("printing peaks file in BED format for ",fileName_chip,"\n")

if(is.null(opt$out_file)){
  fileName <- paste0(opt$output_dir,fileName_chip,".bed")
  write.table(bed,file=fileName,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

}else{
  write.table(bed,file=paste0(opt$output_dir,opt$out_file,".bed"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

}

cat(fileName_chip,"DONE!","\n")
