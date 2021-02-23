args=(commandArgs(TRUE))

genecode<-args[1]      #GTF
bed_ZCL<-args[2]   #bed ZCL
output_NAME<-args[3]     #output filename




library(ChIPpeakAnno, quietly=T)


parseENCODE <- function(x) {

	tmp <- unlist(strsplit(unlist(strsplit(x, ";")), " "))
	tmp <- tmp[tmp != ""]
	tmp <- tmp[seq(2, 12, 2)]
	return(tmp)
}

### AHORA CARGAMOS GENCODE y SELECCIONAMOS LAS COLUMNAS QUE NOS INTERESAN PARA ESTE TIPO DE ANALISIS


genecodev25_all <- read.table(genecode, skip = 5, header = FALSE, sep = "\t")

genecodev25 <- genecodev25_all[genecodev25_all$V3 == "gene", ]
genecodev25_tmp <- apply(as.data.frame(genecodev25[,9]), 1, parseENCODE)
genecodev25_tmp <- t(genecodev25_tmp)

colnames(genecodev25_tmp) <- c("gene_id", "gene_type", "gene_name", "gene_status", "level", "havana_gene")
genecodev25_Annot <- data.frame(genecodev25[,1:8], genecodev25_tmp[, c(1, 2,3, 4, 5,6)])

genecodev25_Annot_Gene <- genecodev25_Annot

G25AnnotChIP <- unique(genecodev25_Annot[genecodev25_Annot[, "V3"] == "gene", c("V1", "V4", "V5", "V7", "gene_id")])
colnames(G25AnnotChIP) <- c("chr", "start", "end", "strand", "gene_id")
G25AnnotChIP_RL <- RangedData(IRanges(start = G25AnnotChIP$start, end = G25AnnotChIP$end, names = paste(G25AnnotChIP$gene_id)), space = paste(G25AnnotChIP$chr), strand = paste(G25AnnotChIP$strand))

#YA TENEMOS EL OBJETO genecodev25_Annot_Gene que contiene la información que nos interesa de Gencode25 y el objeto G25AnnotChIP_RL que es una estructura de datos que contiene a qué gen corresponde cada rango de secuencias. Necesitamos los dos

#COPIA Y PEGA LA FUNCION annotChIP

annotChIP<- function(file_bed, file_annot, file_annot_gene) {

	beddata_RL <- BED2RangedData(file_bed, header=FALSE)
	beddata_Annot <- annotatePeakInBatch(beddata_RL, AnnotationData =file_annot, output = "both", multiple = F, maxgap = 0)
	beddata_Annot.df <- as.data.frame(beddata_Annot)
	beddata_Annot_G15.df <- merge(beddata_Annot.df, file_annot_gene, by.x = 7, by.y = 9)
	beddata.Filter10 <- beddata_Annot_G15.df[((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "upstream") & (abs(beddata_Annot_G15.df$distancetoFeature)<10000)) | ((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "inside") & (abs(beddata_Annot_G15.df$distancetoFeature)<10000)) | ((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "overlapStart")) | ((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "includeFeature")),]

	return(beddata.Filter10)
}

#ahora cargamos el fichero de salida de MACS2 de la actividad 2

data_ZCL<-read.table(bed_ZCL)

#eliminamos las columnas no necesarias y nos quedamos con las que conformarían un fichero *BED

ZCL_bed<-data_ZCL[,c("V1","V2","V3")]

#llamamos a la función annotChIP

bed_data_annotated<-annotChIP(ZCL_bed,G25AnnotChIP_RL,genecodev25_Annot_Gene)

#### Escribimos la salida a un fichero de texto
write.table(bed_data_annotated, output_NAME, quote = FALSE, col.names=T, row.names = FALSE, sep="\t")

#fin del ejercicio
