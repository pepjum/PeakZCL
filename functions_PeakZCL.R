require(wmtsa, quietly = TRUE)
require(multtest,quietly = TRUE)
require(wavelets,quietly = TRUE)
require(Rwave,quietly = TRUE)
require(wavethresh,quietly = TRUE)
require(fields,quietly = TRUE)
require(signal,quietly = TRUE)
require(data.table,quietly = TRUE)
require(VennDiagram,quietly = TRUE)
require(optparse,quietly = TRUE)
require(IRanges,quietly = TRUE)


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

	CI <- findzerocrossing((-1)*amplitude)

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

			if ((poscurr>5) & (poscurr<(s2-5)) & (i>5)) {
				if (length(sum(CI[i, (poscurr-5):(poscurr+5)])) > 0) {
					c7 <- sum(CI[i, (poscurr-5):(poscurr+5)])
				} else {
					c7 <- 0
				}
			} else {
				c7 <- 0
			}

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

			} else if ((poscurr>5) & (poscurr<(s2-5)) & (i>5) & (c7 == 1)) {

					poscurr <- poscurr + sum(mvec11 * CI[i, (poscurr-5):(poscurr+5)])
					klen[kl+1] <- poscurr
					kl <- kl+1
					CI[i, poscurr] <- 0

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
		fc <-sum(sample)/sum(control)
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

peak_calling_ZCL<-function(ZCL,decimate_signal,levth,threshold){

	if (levth == 0 | max(unlist(ZCL$zeroline)) < levth) {
		hLevel <- unlist(ZCL$zeroline)
		hLevel <- hLevel[hLevel>0]
		levth <- (1-cumsum(table(hLevel)/sum(table(hLevel))))/max(1-cumsum(table(hLevel)/sum(table(hLevel))))<0.3
		levth <- which(levth)[1]
		levth<-as.numeric(levth)
		candidatePeaksPos <- which((unlist(ZCL$zeroline) > levth ))

	}
	else{
		candidatePeaksPos <- which((unlist(ZCL$zeroline) > levth))
	}

	return(candidatePeaksPos)
}

peaks_start_end_ZCL<-function(candidatePeaksPos,decimate_signal,factor_decimation){
	putPeaks<- vector("numeric", length(decimate_signal))
	putPeaks[candidatePeaksPos]<-1   #en las posiciones donde están los picos mete un 1
	positions_peaks_center<-which(putPeaks ==1)
	positions_peaks_estimated_left<-sapply(positions_peaks_center,FUN=function(x) x -(5000/factor_decimation))
	positions_peaks_estimated_right<-sapply(positions_peaks_center,FUN=function(x) x +(5000/factor_decimation))

	noisy_signal<-decimate_signal
	for(i in 1:length(positions_peaks_estimated_left)){
		noisy_signal[positions_peaks_estimated_left[i]:positions_peaks_estimated_right[i]]<-999999999999

	}
	selected_signal<-noisy_signal[-(which(noisy_signal ==999999999999))]

	putPeaks[(which(putPeaks == 0 ))] <- 2   #cambia los 0 por 2 no se para que. Lo sigo
	lessthanth <- which(decimate_signal <= mean(selected_signal))
	#lessthanth<-which(decimate_signal <= threshold)
	putPeaks[lessthanth] <- 0                    # menor que threshold a 0
	putPeaks.left <- c(0, putPeaks[-length(putPeaks)])
	putPeaks.right <- c(putPeaks[-1], 0)

	tmp <- putPeaks
	indTmp <- which((putPeaks.left == 2) & (putPeaks == 1) & (putPeaks.right == 2))
	tmp[indTmp-1] <- 1
	tmp[indTmp+1] <- 1

	npeaksTmpPre <- sum(tmp == 2)
	npeaksTmpPost <- 0
		print(paste0(npeaksTmpPre, " Vs ", npeaksTmpPost))

	while (npeaksTmpPre != npeaksTmpPost) {
		npeaksTmpPre <- npeaksTmpPost
		putPeaks.left <- c(0, tmp[-length(tmp)])
		putPeaks.right <- c(tmp[-1], 0)
		indTmpL <- which((putPeaks.left == 2) & (tmp == 1))
		indTmpR <- which((tmp == 1) & (putPeaks.right == 2))
		tmp[indTmpL-1] <- 1
		tmp[indTmpR+1] <- 1
		npeaksTmpPost <- sum(tmp == 2)
	}

	tmp[tmp == 2] <- 0

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


clustering<-function(dataframe,clustering_limit){
	list_sorted<-sort(c(dataframe$start,dataframe$end))      #length
	putative_peaks<- data.frame("start" = list_sorted[seq(1, nchar(list_sorted), 2)], "end" = list_sorted[seq(2, nchar(list_sorted), 2)])

	join<-F
  new <- 0
	old<-nrow(dataframe)
  while ((new) != (old)){
		list_sorted<-sort(c(putative_peaks$start,putative_peaks$end))
		putative_peaks<- data.frame("start" = list_sorted[seq(1, nchar(list_sorted), 2)], "end" = list_sorted[seq(2, nchar(list_sorted), 2)])
		start_peaks_clustered<-data.frame("start"=NULL,"end"=NULL)
		print(paste0(new, "  Vs  ", old))

  	for (i in 1:(nrow(putative_peaks)-1)){
  			if(join){
  				join<-F
  				next
  				}else{
  					if((putative_peaks$start[i+1])-(putative_peaks$end[i]) <= clustering_limit ){

  						start_peaks_clustered <- rbind(start_peaks_clustered, data.frame(start=putative_peaks$start[i],end=putative_peaks$end[i+1]))

  						join<-T
  						}else{
							start_peaks_clustered <- rbind(start_peaks_clustered, data.frame(start=putative_peaks$start[i],end=putative_peaks$end[i]))

  						}

  					}

  				}
					n_putative_peaks<-nrow(putative_peaks)
					number<-nrow(start_peaks_clustered)
					print(number)
					if((putative_peaks$start[n_putative_peaks])-(start_peaks_clustered$end[number]) <= clustering_limit )
					{
						start_peaks_clustered<-rbind(start_peaks_clustered,data.frame(start=start_peaks_clustered$start[number],end=putative_peaks$end[n_putative_peaks]))
					}else{
						start_peaks_clustered<-rbind(start_peaks_clustered,data.frame(start=putative_peaks$start[n_putative_peaks],end=putative_peaks$end[n_putative_peaks]))
					}

					old <- nrow(putative_peaks)
					new <- nrow(start_peaks_clustered)
					putative_peaks<-start_peaks_clustered
        }

				return(putative_peaks)
}

bining<-function(dataframe,Dalton){
    list_sorted<-sort(c(dataframe$start,dataframe$end))
    dataframe<- data.frame("M" = list_sorted[seq(1, length(list_sorted), 2)], "Z" = list_sorted[seq(2, length(list_sorted), 2)])

    #dataframe<-as.data.frame(dataframe)
    #dataframe<-dataframe[!(dataframe$M=="NA"),]
	new = 0
	old = nrow(dataframe)

	dataframe <- setNames(dataframe,c("M","Z"))
	while( old != new){
        print(paste0(new, "  Vs  ", old))
        list_sorted<-sort(c(dataframe$M,dataframe$Z))
        dataframe<- data.frame("M" = list_sorted[seq(1, length(list_sorted), 2)], "Z" = list_sorted[seq(2, length(list_sorted), 2)])
		tmp = data.frame('M' = NULL, 'Z' = NULL)
		SKIP_LINE = F
		iteration = nrow(dataframe)- 1
		for (i in 1:(iteration)){
			if(SKIP_LINE){
				SKIP_LINE = F
				next
			}
			dataframe<-dataframe[order(as.numeric(paste(dataframe$M))),]
			if(as.numeric(paste(dataframe[i+1, 'M'])) - as.numeric(paste(dataframe[i, 'Z']))  <= as.numeric(Dalton)){
				#print('shrinking')
				newM = as.numeric(paste(dataframe[i, 'M']))
				newZ = as.numeric(paste(dataframe[i+1, 'Z']))
				tmp = rbind(tmp, data.frame('M' = as.numeric(newM), 'Z' = as.numeric(newZ)))
				SKIP_LINE = T
			}else{
				tmp = rbind(tmp, dataframe[i,])
				SKIP_LINE = F
			}
		}
		tmp <- rbind(tmp, dataframe[i+1,])
		new <- nrow(tmp)
		old <- nrow(dataframe)
		dataframe <- tmp
	}
	if(as.numeric(paste(dataframe[nrow(dataframe), 'M'])) - as.numeric(paste(dataframe[nrow(dataframe)-1, 'Z'])) <= as.numeric(Dalton)){
        newM = as.numeric(paste(dataframe[i, 'M']))
        newZ = as.numeric(paste(dataframe[i+1, 'Z']))
		dataframe2 <- dataframe[c(seq(1,nrow(dataframe)-2)), ]
		dataframe2 <- rbind(dataframe2, data.frame('M' = as.numeric(newM), 'Z' = as.numeric(newZ)))
	}else{
        dataframe2 <- dataframe

    }
    dataframe2<-unique(dataframe2)
	names(dataframe2)<-c("start","end")
	return(dataframe2)

}

quantification<-function(dataframe,signal){
  quantified_peaks<-list()
  for (i in 1:nrow(dataframe)){
    peaks_quantified<-sum(signal[dataframe$start[i]:dataframe$end[i]])
    quantified_peaks[[i]]<-peaks_quantified
  }
  return(quantified_peaks)
}


plot_ZCL_lineTRACE<-function(ZCL){
	cat("generating plot ZCL....","\n")
	for (i in 1:length(ZCL$zerolineTrace)){
		if (unlist(ZCL$zeroline)[i]>0) {
			x<-c(unlist(ZCL$zerolineTrace[[i]]))
			y<-c(seq(1:length(x)))
			if (i == which(unlist(ZCL$zeroline) != 0)[1]) {
				plot(x,y, type="l", xlim = c(1, 5001), ylim = c(1, 50))
				} else {
					lines(x,y, type="l")
				}
			}
		}
		dev.off()
}

estimating_CWT_and_peak_calling<-function(signal,length_vector_wavelet,scale_wavelets,levth_number,factor_decimation){

    if(length(signal)>length_vector_wavelet){
      list_wavelets<-split(signal, ceiling(seq_along(signal)/length_vector_wavelet))
			#cat(length(list_wavelets),"- parts splitted signal","\n")
		peaks_chip_list<-list()
      for(i in 1:length(list_wavelets)){
				#signal<-unlist(list_wavelets[i])
        CWT_signal<-wavCWT(list_wavelets[[i]], n.scale=scale_wavelets,wavelet="gaussian1")
				#cat("Zero-crossing lines",i,"\n")
        ZCL<-zeroCrossingLineTrace(t(as.matrix(CWT_signal)))
        #cat("peak calling",i,"\n")
				peaks_chip_list[[i]]<-peak_calling_ZCL(ZCL,list_wavelets[[i]],levth=levth_number)
				peaks_chip_list[[i]] <- ((i-1)*length_vector_wavelet)+peaks_chip_list[[i]]
			}
			peaks_chip<-unlist(peaks_chip_list)

      #cat("finding start and end of peaks","\n")
      start_end_peaks<-peaks_start_end_ZCL(peaks_chip,signal,factor_decimation)
    }else{
			CWT_signal<-wavCWT(signal, n.scale=scale_wavelets,wavelet="gaussian1")
      #cat("Zero-crossing lines...","\n")
      ZCL<-zeroCrossingLineTrace(t(as.matrix(CWT_signal)))
      #cat("peak calling....","\n")
			peaks_chip<-peak_calling_ZCL(ZCL,signal,levth=levth_number)
			start_end_peaks<-peaks_start_end_ZCL(peaks_chip,signal,factor_decimation)
  	}
  return(start_end_peaks)
}

statistics<-function(peaks_in_position,area_chip_number,fold_change_number,quantified,quantified_input){
  peaks_in_position$area_chip<-log10(unlist(quantified))
  peaks_in_position$area_input<-log10(unlist(quantified_input))
  peaks_in_position$fold_change<-(peaks_in_position$area_chip)-(peaks_in_position$area_input)
  nlocation_new<-1:length(peaks_in_position$fold_change)
  n_location_inf<-nlocation_new[peaks_in_position$fold_change == Inf]
  peaks_in_position$fold_change[n_location_inf]<-peaks_in_position$area_chip[n_location_inf]
  peaks_in_position$zscore<-scale(peaks_in_position$fold_change)
  #peaks_in_position$abundancia<-((log10(unlist(quantified)))+(log10(unlist(quantified_input)))/2)
	peaks_in_position$selected <- ((peaks_in_position$area_chip > area_chip_number) & (peaks_in_position$fold_change > fold_change_number))*1
  return(peaks_in_position)#peaks_results<-peaks_in_position[which(peaks_in_position$zscore >=1.64),]
}

BED_files<-function(dataframe,chr_number,filePath_control,all_peaks){

	if(!is.null(filePath_control)){
			df_selected<-dataframe[which(dataframe$selected ==1),]
			score<-df_selected$zscore
			chr<-seq(1:nrow(df_selected))
			pico<-seq(1:nrow(df_selected))
			for (i in 1:length(chr)){
	    	chr[i]<-chr_number
				pico[i]<-i
	  	}
  		bed<-data.frame("chr"=as.character(unlist(paste0(chr))),"start"=as.numeric(df_selected$start),"end"=as.numeric(df_selected$end),"name"=as.character(paste0("A_Chip:",as.numeric(df_selected$area_chip)," A_Input:",as.numeric(df_selected$area_input)," fchange:",as.numeric(df_selected$fold_change)),"score"=as.numeric(unlist(score))))

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
