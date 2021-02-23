### model

args=(commandArgs(TRUE))


Bed<-args[1]      #Bed file normalization
filter<-args[2]
peaks_clustered<-read.table(Bed)
names(peaks_clustered)<-c("chr","start","end","name","counts_chip","counts_input")
#names(data)<-c("chr","start","end","name","counts_chip","counts_input","area_chip","area_input","fold_change_log2")

# fold_change_order<-sort(peaks_clustered$fold_change_log2)
#
# file_new_reorder<-peaks_clustered[order(peaks_clustered$fold_change_log2, decreasing=T),]
#
#
# fold_change_order_finite<-fold_change_order[is.finite(fold_change_order)]
#
# df<-data.frame("pico"=seq(1:length(fold_change_order_finite)), "fold_change"=fold_change_order_finite)
#
# point1<-df[which.min(abs(df$fold_change -median(df$fold_change ))),]
#
# point2<-df[which.min(abs(df$fold_change -median(df$fold_change )+0.5*(mad(df$fold_change)))),]
# point3<-df[which.min(abs(df$fold_change -median(df$fold_change )+1*(mad(df$fold_change)))),]
# point4<-df[which.min(abs(df$fold_change -median(df$fold_change )-0.5*(mad(df$fold_change)))),]
# point5<-df[which.min(abs(df$fold_change -median(df$fold_change )-1*(mad(df$fold_change)))),]
#
# recta<-rbind(point3,point2, point1, point4, point5)
#
# reg<-lm(recta$pico ~ recta$fold_change)
#
# coef1<- as.numeric(paste(reg$coefficients[1]))
# coef2<- as.numeric(paste(reg$coefficients[2]))
#
# file_new_reorder$pico_real<-nrow(file_new_reorder):1
# file_new_reorder$pico_estimado<-file_new_reorder$fold_change_log2 * coef2 + coef1
#
# file_new_reorder$diff<-file_new_reorder$pico_estimado - file_new_reorder$pico_real
#
#
# #selected<-file_new_reorder[which(file_new_reorder$diff >  summary(reg)[[6]])]
#
# ########
#
# file_new_reordered<-file_new_reorder[order(file_new_reorder$diff, decreasing=T),]
# rownames(file_new_reordered)<-seq(1:nrow(file_new_reordered))
#
# cut_line<-file_new_reordered[which.min(abs(file_new_reordered$diff - summary(reg)[[6]])),]
#
# selected<-file_new_reordered[1:(rownames(cut_line)),]
#
# selected$areas<-NULL
# #selected$fold_change<-NULL
# selected$diff<-NULL
# selected$pico_real<-NULL
# selected$pico_estimado<-NULL


########

quita0<-peaks_clustered[((peaks_clustered$counts_chip - peaks_clustered$counts_input) >0),]

if (filter=="mean"){
    selected<-quita0[which(quita0$counts_chip > as.numeric(paste(summary(quita0$counts_chip)[4]))),]
}else if (filter=="median"){#selected<-quita0[which(quita0$counts_chip > mean(quita0$counts_chip)),]
    selected<-quita0[which(quita0$counts_chip > as.numeric(paste(summary(quita0$counts_chip)[3]))),]
}


write.table(selected, file=paste0(lapply(strsplit(paste(Bed),".Norm"),"[",1),".model"),col.names=F, row.names=F, quote=F, sep="\t")
