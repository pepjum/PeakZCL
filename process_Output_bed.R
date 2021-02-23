args=(commandArgs(TRUE))


file <- args[1]

library(stringr)

file_l<-read.table(file)

file_l_chrX<-file_l[which(file_l$V1=="chrX"),]
file_l_chrY<-file_l[which(file_l$V1=="chrY"),]
file_l_chrM<-file_l[which(file_l$V1=="chrM"),]

b<-file_l[!(file_l$V1 %in% c("chrX","chrY","chrM")),]

b$number<-str_extract(paste(b$V1),"[0-9]+")
b_ordered <- b[order(as.numeric(b$number), as.numeric(b$V2)),]

b_ordered$number<-NULL

if((nrow(file_l_chrX)>0)){
    b_ordered<-rbind(b_ordered,file_l_chrX)
}
if((nrow(file_l_chrY)>0)){
    b_ordered<-rbind(b_ordered,file_l_chrY)
}
if((nrow(file_l_chrM)>0)){
    b_ordered<-rbind(b_ordered,file_l_chrM)
}



path_output<-dirname(file)

write.table(b_ordered,file=paste0(path_output,"/ZCL_peaks_ordered.bed"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
