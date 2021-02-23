
suppressPackageStartupMessages({
	require(wmtsa)
	require(multtest)
	require(wavelets)
	require(Rwave)
	require(wavethresh)
	require(fields)
	require(signal)
	require(data.table)
	require(VennDiagram)
	require(optparse)
	require(IRanges)

})




parseENCODE <- function(x) {

	tmp <- unlist(strsplit(unlist(strsplit(x, ";")), " "))
	tmp <- tmp[tmp != ""]
	tmp <- tmp[seq(2, 12, 2)]
	return(tmp)
}



annotChIP<- function(file_bed, file_annot, file_annot_gene) {

	beddata_RL <- BED2RangedData(file_bed, header=FALSE)

	beddata_Annot <- annotatePeakInBatch(beddata_RL, AnnotationData =file_annot, output = "both", multiple = F, maxgap = 0)
	beddata_Annot.df <- as.data.frame(beddata_Annot)
	#write.table(beddata_Annot.df, file = paste(nameOut, "_Annot.txt", sep = ""), quote = FALSE, row.names = FALSE, sep="\t")
	beddata_Annot_G15.df <- merge(beddata_Annot.df, file_annot_gene, by.x = 7, by.y = 9)
	#write.table(beddata_Annot_G15.df, file = paste(nameOut, "_Annot_G15.txt", sep = ""), quote = FALSE, row.names = FALSE, sep="\t")

	beddata.Filter10 <- beddata_Annot_G15.df[((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "upstream") & (abs(beddata_Annot_G15.df$distancetoFeature)<10000)) | ((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "inside") & (abs(beddata_Annot_G15.df$distancetoFeature)<10000)) | ((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "overlapStart")) | ((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "includeFeature")),]

	#write.table(beddata.Filter10, file = paste(nameOut, "_Filter10.txt", sep = ""), quote = FALSE, row.names = FALSE, sep="\t")
	return(beddata.Filter10)
}


 findzerocrossing<- function(I) {

	s1 <- nrow(I)
	s2 <- ncol(I)

	crossI <- matrix(0, s1, s2)

	for (i in 1:s1) {
		for (k in 1:(s2-1)) {

			if ((I[i,k]*I[i,k+1])<0) {

				if (I[i,k]<0) {

					crossI[i,k] = 1;

				} else {

					crossI[i,k] = (-1);

				}

			} else if ((I[i,k]==0)&(k>1)) {

				if ((I[i,k-1]<0)&(I[i,k+1]>0)) {

					crossI[i,k] = 1;

				} else if ((I[i,k-1]>0)&(I[i,k+1]<0)) {

					crossI[i,k] = (-1);

				}

			}
		}
	}

	return(crossI)

}

#function that calculates the joints by the cut points between the scales


zeroCrossingLineTrace<- function(amplitude) {

	CI <- abs(findzerocrossing(amplitude))

	s1 <- nrow(amplitude)
	s2 <- ncol(amplitude)
	zerolineTrace <- list()
	length(zerolineTrace) <- s2
	zeroline <- matrix(0, 1, s2)

	mvec11 <- c(-5:5)
	mvec9 <- c(-4:4)
	mvec7 <- c(-3:3)
	mvec5 <- c(-2:2)
	mvec3 <- c(-1:1)

	for (pos in 1:(s2-1)) {

		klen <- list()
		kl <- 0
		poscurr <- pos
		permit <- 0

		for (i in 1:s1) {

			if (length(sum(CI[i, (poscurr-1):(poscurr+1)])) > 0) {
				c1 <- sum(CI[i, (poscurr-1):(poscurr+1)])
			} else {
				c1 <- 0
			}

			if (length(CI[i, poscurr+1]) > 0) {
				c2 <- CI[i, poscurr+1]
			} else {
				c2 <- 0
			}

			if (length(CI[i, poscurr-1]) > 0) {
				c3 <- CI[i, poscurr-1]
			} else {
				c3 <- 0
			}

			if ((poscurr>2) & (poscurr<(s2-2)) & (i>2)) {
				if (length(sum(CI[i, (poscurr-2):(poscurr+2)])) > 0) {
					c4 <- sum(CI[i, (poscurr-2):(poscurr+2)])
				} else {
					c4 <- 0
				}
			} else {
				c4 <- 0
			}

			if ((poscurr>3) & (poscurr<(s2-3)) & (i>3)) {
				if (length(sum(CI[i, (poscurr-3):(poscurr+3)])) > 0) {
					c5 <- sum(CI[i, (poscurr-3):(poscurr+3)])
				} else {
					c5 <- 0
				}
			} else {
				c5 <- 0
			}

			if ((poscurr>4) & (poscurr<(s2-4)) & (i>4)) {
				if (length(sum(CI[i, (poscurr-4):(poscurr+4)])) > 0) {
					c6 <- sum(CI[i, (poscurr-4):(poscurr+4)])
				} else {
					c6 <- 0
				}
			} else {
				c6 <- 0
			}

			# if ((poscurr>5) & (poscurr<(s2-5)) & (i>5)) {
			# 	if (length(sum(CI[i, (poscurr-5):(poscurr+5)])) > 0) {
			# 		c7 <- sum(CI[i, (poscurr-5):(poscurr+5)])
			# 	} else {
			# 		c7 <- 0
			# 	}
			# } else {
			# 	c7 <- 0
			# }

			if (CI[i, poscurr] == 1) {

				klen[kl+1] <- poscurr
				kl <- kl+1
				CI[i, poscurr] <- 0

			} else if ((poscurr>1) & (poscurr<s2-1) & (i>1) & (c1 == 1)) {

					poscurr <- poscurr + sum(mvec3 * CI[i, (poscurr-1):(poscurr+1)])
					klen[kl+1] <- poscurr
					kl <- kl+1
					CI[i, poscurr] <- 0

			} else if ((poscurr == 1) & (i>1) & (c2 == 1)) {

				poscurr <- poscurr+1
				klen[kl+1] <- poscurr
				kl <- kl+1
				CI[i, poscurr] <- 0

			} else if ((poscurr == (s2-1)) & (i>1) & (c3 == 1)) {

					poscurr <- poscurr-1
					klen[kl+1] <- poscurr
					kl <- kl+1
					CI[i, poscurr] <- 0

			} else if ((poscurr>2) & (poscurr<(s2-2)) & (i>2) & (c4 == 1)) {

					poscurr <- poscurr + sum(mvec5 * CI[i, (poscurr-2):(poscurr+2)])
					klen[kl+1] <- poscurr
					kl <- kl+1
					CI[i, poscurr] <- 0

			} else if ((poscurr>3) & (poscurr<(s2-3)) & (i>3) & (c5 == 1)) {

					poscurr <- poscurr + sum(mvec7 * CI[i, (poscurr-3):(poscurr+3)])
					klen[kl+1] <- poscurr
					kl <- kl+1
					CI[i, poscurr] <- 0

			} else if ((poscurr>4) & (poscurr<(s2-4)) & (i>4) & (c6 == 1)) {

					poscurr <- poscurr + sum(mvec9 * CI[i, (poscurr-4):(poscurr+4)])
					klen[kl+1] <- poscurr
					kl <- kl+1
					CI[i, poscurr] <- 0

			# } else if ((poscurr>5) & (poscurr<(s2-5)) & (i>5) & (c7 == 1)) {
            #
			# 		poscurr <- poscurr + sum(mvec11 * CI[i, (poscurr-5):(poscurr+5)])
			# 		klen[kl+1] <- poscurr
			# 		kl <- kl+1
			# 		CI[i, poscurr] <- 0

			} else {

				if ((permit>1) | (i == 1)) break;
				permit <- permit+1

			}

		}

		if (length(klen) > length(zerolineTrace[[pos]])) {

			zerolineTrace[pos] <- list(unlist(klen))
			zeroline[pos] <- list(kl)

		}

	}

	return(list(zerolineTrace = zerolineTrace, zeroline = zeroline))

}

NormalizationChIP<- function( sample, control, method){

	if (method == "SDS"){
		fc <-sum(as.numeric(sample))/sum(as.numeric(control))
	}
	else if (method == "SES"){
		Y <- sample
		X <- control
		Ysorted<- sort(Y)  #sort ordena de menor a may
		Xsorted <- X[order(Y)]
		k <- which.max(abs((cumsum(Ysorted)/sum(Y))-(cumsum(Xsorted)/sum(X)))) #k es una posicion
		Yk <- cumsum(Ysorted)[k]
		Xk <- cumsum(Xsorted)[k]
		fc <- Yk/Xk
	}
	return(fc)
}
# outputs<-data.frame("Ysorted"=Ysorted)
# cumsums<-apply(outputs, 2, cumsum)
#
# plot(Yk)

peak_calling_ZCL<-function(ZCL,decimate_signal,levth, noiselevel){

	if (levth == 0 | max(unlist(ZCL$zeroline)) < levth) {
		hLevel <- unlist(ZCL$zeroline)
		hLevel <- hLevel[hLevel>0]
		levth <- (1-cumsum(table(hLevel)/sum(table(hLevel))))/max(1-cumsum(table(hLevel)/sum(table(hLevel))))<0.3
		levth <- which(levth)[1]
		levth<-as.numeric(levth)

		#candidatePeaksPos_sin<-which((unlist(ZCL$zeroline) > 0))
		peaks<-unlist(ZCL$zeroline)
		peaks[decimate_signal < median(decimate_signal)] = 0
		candidatePeaksPos <- which(peaks > levth)

	}
	else{
		#candidatePeaksPos_sin<-which((unlist(ZCL$zeroline) >= 0))
		peaks<-unlist(ZCL$zeroline)
		peaks[decimate_signal < median(decimate_signal)] = 0
		candidatePeaksPos <- which(peaks > as.numeric(paste(levth)))
	}

	intensities <- decimate_signal[unlist(lapply(ZCL$zerolineTrace[candidatePeaksPos], FUN = function(x) x[1]))]

#	peaks_candidates<-candidatePeaksPos[decimate_signal[unlist(lapply(ZCL$zerolineTrace[(candidatePeaksPos)], FUN = function(x) x[1]))] > median(intensities)]
	peaks_candidates<-candidatePeaksPos[decimate_signal[unlist(lapply(ZCL$zerolineTrace[(candidatePeaksPos)], FUN = function(x) x[1]))] > as.numeric(paste(noiselevel))]

	#peaks_candidates_nofiltered_longitudes<-candidatePeaksPos_sin[decimate_signal[unlist(lapply(ZCL$zerolineTrace[(candidatePeaksPos_sin)], FUN = function(x) x[1]))]>0]
	#peaks_candidates_nofiltered_intensities<-candidatePeaksPos[decimate_signal[unlist(lapply(ZCL$zerolineTrace[(candidatePeaksPos)], FUN = function(x) x[1]))] > 0]

	return(peaks_candidates)
	#return(peaks_candidates_nofiltered_intensities)
}

peaks_start_end_ZCL<-function(candidatePeaksPos,decimate_signal,threshold, noiselevel){
	putPeaks<- vector("numeric", length(decimate_signal))
	putPeaks[candidatePeaksPos]<-1   #en las posiciones donde están los picos mete un 1
	putPeaks[(which(putPeaks == 0 ))] <- 2   #cambia los 0 por 2 no se para que. Lo sigo
	#lessthanth <- which(decimate_signal <= summary(decimate_signal)[[threshold]])
	lessthanth <- which(decimate_signal <= noiselevel)

	#lessthanth<-which(decimate_signal <= threshold)
	putPeaks[lessthanth] <- 50                    # menor que threshold a 0
	#putPeaks[lessthanth] <-2
	putPeaks.left <- c(0, putPeaks[-length(putPeaks)]) # le mete un 0 delante y otro detras para corregir
	putPeaks.right <- c(putPeaks[-1], 0)

	tmp <- putPeaks
	indTmp <- which((putPeaks.left == 2 ) & (putPeaks == 1) & (putPeaks.right == 2 ))
	#indTmp <- which((putPeaks.left == 2 || putPeaks.left == 0 ) & (putPeaks == 1) & (putPeaks.right == 2 || putPeaks.right==0))

	tmp[indTmp-1] <- 1
	tmp[indTmp+1] <- 1

	npeaksTmpPre <- sum(tmp == 2)
	npeaksTmpPost <- 0
		#print(paste0(npeaksTmpPre, " Vs ", npeaksTmpPost))

	while (npeaksTmpPre != npeaksTmpPost) {
		npeaksTmpPre <- npeaksTmpPost
		putPeaks.left <- c(0, tmp[-length(tmp)])
		putPeaks.right <- c(tmp[-1], 0)
		#indTmpL <- which((putPeaks.left == 2 ) & (tmp == 1))
		indTmpL <- which((putPeaks.left == 2 ) & (tmp == 1))

		#indTmpR <- which((tmp == 1) & (putPeaks.right == 2 ))
		indTmpR <- which((tmp == 1) & (putPeaks.right == 2 ))

		tmp[indTmpL-1] <- 1
		tmp[indTmpR+1] <- 1
		npeaksTmpPost <- sum(tmp == 2)
	}

	tmp[tmp == 2] <- 0
	tmp[tmp == 50] <-0

	if(tmp[1] == 1){
		tmp[1] <- 0
		tmp[length(tmp)] <- 0
	}

	peaksPos <- tmp                         #yo ya tengo los picos en su posicion, no tengo q recuperar 0

	start_peak <- which(diff(peaksPos) == 1)+1
	end_peak <- which(diff(peaksPos) == (-1))+1

	putative_peaks <- data.frame("start" = start_peak, "end"= end_peak) #"summit"=candidatePeaksPos)
	return(putative_peaks)
}

quantification<-function(dataframe,signal){
  quantified_peaks<-list()
  for (i in 1:nrow(dataframe)){
    peaks_quantified<-sum(signal[dataframe$start[i]:dataframe$end[i]])
    quantified_peaks[[i]]<-peaks_quantified
  }
  return(quantified_peaks)
}


plot_ZCL_lineTRACE<-function(ZCL, file_output, scales){
	cat("generating plot ZCL....","\n")
	pdf(file_output,,width=20, height=14, onefile=FALSE)
	for (i in 1:length(ZCL$zerolineTrace)){
		if (unlist(ZCL$zeroline)[i]>0) {
			x<-c(unlist(ZCL$zerolineTrace[[i]]))
			y<-c(seq(1:length(x)))
			if (i == which(unlist(ZCL$zeroline) != 0)[1]) {
				plot(x,y, type="l", xlim = c(1, length(ZCL$zerolineTrace)), ylim = c(1, scales))
				} else {
					lines(x,y, type="l")
				}
			}
		}
	dev.off()
}

estimating_CWT_and_peak_calling<-function(signal,length_vector_wavelet,scale_wavelets,levth_number,decimating,threshold, noiselevel){

    if(length(signal)>length_vector_wavelet){
      list_wavelets<-split(signal, ceiling(seq_along(signal)/length_vector_wavelet))

		peaks_chip_list<-list()
      for(i in 1:length(list_wavelets)){

        CWT_signal<-wavCWT(list_wavelets[[i]],scale.range=c(1,scale_wavelets), n.scale=scale_wavelets,wavelet="gaussian1")
				#cat("Zero-crossing lines",i,"\n")
        ZCL<-zeroCrossingLineTrace(t(as.matrix(CWT_signal)))
        #cat("peak calling",i,"\n")
				peaks_chip_list[[i]]<-peak_calling_ZCL(ZCL,list_wavelets[[i]],levth=levth_number, noiselevel)
				peaks_chip_list[[i]] <- ((i-1)*length_vector_wavelet)+peaks_chip_list[[i]]
		}
		peaks_chip<-unlist(peaks_chip_list)

      #cat("finding start and end of peaks","\n")
      start_end_peaks<-peaks_start_end_ZCL(peaks_chip,signal,threshold, noiselevel)
    }else{
			CWT_signal<-wavCWT(signal, scale.range=c(1,scale_wavelets), n.scale=scale_wavelets,wavelet="gaussian1")
      #cat("Zero-crossing lines...","\n")
      		ZCL<-zeroCrossingLineTrace(t(as.matrix(CWT_signal)))
      #cat("peak calling....","\n")
			peaks_chip<-peak_calling_ZCL(ZCL,signal,levth=levth_number, noiselevel)
			start_end_peaks<-peaks_start_end_ZCL(peaks_chip,signal,threshold, noiselevel)
  	}
  return(start_end_peaks)
}

statistics<-function(peaks_in_position,area_chip_number,fold_change_number,quantified,quantified_input, enrichment){
  peaks_in_position$quant_chip_N<-unlist(quantified)
  peaks_in_position$quant_input_N<-unlist(as.numeric(paste(quantified_input)))
  #peaks_in_position$quant_chip_Norm<- (as.numeric(paste(peaks_in_position$quant_chip_no_N)) / total_counts_chip)
  #peaks_in_position$quant_input_Norm<- (as.numeric(paste(peaks_in_position$quant_input_no_N)) / total_counts_input)

  peaks_in_position<-peaks_in_position[!(peaks_in_position$quant_chip_N==0),]
  peaks_in_position$area_chip<-log2(peaks_in_position$quant_chip_N)
  peaks_in_position$area_input<-log2(peaks_in_position$quant_input_N)
  peaks_in_position$area_input[!is.finite(peaks_in_position$area_input)]<-0

  peaks_in_position$fold_change<-(peaks_in_position$area_chip)-(peaks_in_position$area_input)
  peaks_in_position$zscore<-scale(peaks_in_position$fold_change)

  #peaks_in_position$zscore<-scale(peaks_in_position$fold_change, center = FALSE, scale = apply(x, 2, sd, na.rm = TRUE))


  peaks_in_position$abundancia<-(peaks_in_position$area_chip + peaks_in_position$area_input)/2

	#plotter<-data.frame("fold_change"=paste(peaks_in_position$fold_change), "abundancia"=paste(peaks_in_position$abundancia))
	#peaks_in_position$selected <- ((peaks_in_position$area_chip > area_chip_number) & (peaks_in_position$fold_change > fold_change_number) & (peaks_in_position$abundancia > mean(as.numeric(paste(peaks_in_position$abundancia))-2*(sd(as.numeric(paste(peaks_in_position$abundancia)))))) & (peaks_in_position$intensity_chip > abs(mean(as.numeric(paste(peaks_in_position$intensity_chip))-2*(sd(as.numeric(paste(peaks_in_position$intensity_chip))))))) & (peaks_in_position$width > abs(mean(as.numeric(paste(peaks_in_position$width)))-2*(sd(as.numeric(paste(peaks_in_position$width)))))))*1
	if (enrichment=="histone"){
		peaks_in_position$selected <- ((peaks_in_position$area_chip > area_chip_number) & (peaks_in_position$fold_change > 0.1 ) & (peaks_in_position$abundancia > abs(mean(as.numeric(paste(peaks_in_position$abundancia))-2*(sd(as.numeric(paste(peaks_in_position$abundancia))))))))*1

	}else if (enrichment =="TF"){
#		peaks_in_position$selected <- ((peaks_in_position$area_chip > area_chip_number) & (peaks_in_position$fold_change > fold_change_number) & (peaks_in_position$abundancia > abs(mean(as.numeric(paste(peaks_in_position$abundancia))-2*(sd(as.numeric(paste(peaks_in_position$abundancia))))))) & (peaks_in_position$width > abs(mean(as.numeric(paste(peaks_in_position$width)))-2*(sd(as.numeric(paste(peaks_in_position$width)))))))*1
		peaks_in_position$selected <- ((peaks_in_position$area_chip > area_chip_number) & (peaks_in_position$fold_change > 0.1) & (peaks_in_position$abundancia > abs(mean(as.numeric(paste(peaks_in_position$abundancia))-2*(sd(as.numeric(paste(peaks_in_position$abundancia))))))) )*1


	}

  return(peaks_in_position)#peaks_results<-peaks_in_position[which(peaks_in_position$zscore >=1.64),]
}

BED_files<-function(dataframe,chr_number,filePath_control,all_peaks){

	if(!is.null(filePath_control)){
		if(all_peaks=="FALSE"){
			df_selected<-dataframe[which(dataframe$selected ==1),]
		}else{
			df_selected<-dataframe
		}# score<-df_selected$zscore
			chr<-seq(1:nrow(df_selected))
			pico<-seq(1:nrow(df_selected))
			for (i in 1:length(chr)){
	    	chr[i]<-chr_number
				pico[i]<-i
	  	}
  		bed<-data.frame("chr"=as.character(unlist(paste0(chr))),"start"=as.numeric(df_selected$start),"end"=as.numeric(df_selected$end),"name"=as.character(paste0("A_Chip:",as.numeric(df_selected$area_chip)," A_Input:",as.numeric(df_selected$area_input)," zscore:",as.numeric(df_selected$zscore))))

		}else{
			chr<-seq(1:nrow(dataframe))
			pico<-seq(1:nrow(dataframe))
	  	for (i in 1:length(chr)){
	    	chr[i]<-chr_number
				pico[i]<-i
	  	}
			bed<-data.frame("chr"=as.character(unlist(paste0(chr))),"start"=as.numeric(dataframe$start),"end"=as.numeric(dataframe$end),"name"=as.character(paste0("A_Chip:",as.numeric(dataframe$area_chip))))
		}
  	return(bed)
}
