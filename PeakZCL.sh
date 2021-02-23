#!/bin/bash

printf "\n\n===================================================================\n PeakZCL. A peak caller tool based on Zero-Crossing Lines for ChIP-seq\n===================================================================\n\n"


day=`date +"%d/%m/%Y"`
hour=`date +"%H:%M"`
printf "Started at $day $hour \n"

sPath="`dirname \"$0\"`"
sPath="`( cd \"$sPath\" && pwd )`"

#printf $sPath\n

usage()
{
cat << EOF

OPTIONS:
   -f   Directory containing Sample file (required)
   -c   Directory containing control file
   -n   Normalization method ("SES","SDS", "N", default: SES)
   -w   Decimating signal factor (default: 50)
   -l   Length signal wavelet (default: 200000)
   -s   Scales of the wavelet (default: 40)
   -g   Zero-crossing lines leverage threshold (default: 6)
   -z   Nature of enrichment. Options are ("histone" or "TF". Default: "TF")
   -j   Bp range between peaks to clustering purposes (default: 50)
   -a   Minimum value of the area of the selected peaks given in log10(x) (default: 1.2)
   -k   Difference between peak area of sample and input. The value is given in log10(x) (default : 0.1)
   -o   Output BED file name (default: ZCL_peaks.final.bed)
   -v   Path to annotation database in .gtf format (default: NULL)
   -m   Returns all peaks detected or not (VALUES: TRUE, FALSE. Default=FALSE)
   -p   Number of processors used by R scripts (default: 1)
EOF
}

#Defaults --

normalization="SES"
decimating="50"
length="200000"
scales="40"
levth="6"
clustering="50"
area="1.2"
foldchange="0.1"
cores="1"
sdir=""
cdir=""
odir="NULL"
outnamefile="PeakZCL_peaks"
annotfile=""
enrichment="TF"
allpeaks="FALSE"
while getopts "f:c:n:w:l:s:g:z:j:a:k:o:v:y:m:p:" OPTION
do
	case $OPTION in
	f) sdir=$OPTARG
	;;
	c) cdir=$OPTARG
	;;
	n) normalization=$OPTARG
	;;
	w) decimating=$OPTARG
	;;
	l) length=$OPTARG
	;;
	s) scales=$OPTARG
	;;
	g) levth=$OPTARG
	;;
	z) enrichment=$OPTARG
	;;
	j) clustering=$OPTARG
	;;
	a) area=$OPTARG
	;;
	k) foldchange=$OPTARG
	;;
	o) outnamefile=$OPTARG
	;;
    m) allpeaks=$OPTARG
    ;;
	v) annotfile=$OPTARG
	;;
  p) cores=$OPTARG
  ;;
	?)
	usage
	exit
	;;
	esac
done

if [[ -z $sdir ]]
then
     usage
     exit 1
fi

echo -e $sdir\n
echo -e $cdir\n
#echo -e $cdir\n
#echo -e $chromsizes\n
echo -e "preprocessing chip and input files\n"

EXP=$sdir'chr_names.txt'

if [[ ! -z $sdir'TMP' ]]; then
    mkdir $sdir'TMP'
    mkdir $sdir'TMP/chip/'
fi

if [[ ! -z $cdir ]]; then
  mkdir $sdir'TMP/input/'

fi
#wdir= $sdir/TMP

#if [[ $odir == "NULL" ]]; then
mkdir $sdir'OUTPUT_6/'
odir=$sdir'OUTPUT_6/'
find $odir -maxdepth 1 -type f -delete

#fi
export LANG=C
export LC_ALL=C

#preprocessing chip
if [[ ! -z $sdir ]]; then
  printf "creating intermediate chip files\n"

  while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line
    chrom=$line
    #printf "$line\n"
    # echo -e  "grep "$chrom\>" $sdir/*.bedGraph >> $sdir/TMP/chip/$line.chrom &"
    grep "$chrom\>" $sdir/*.bedGraph >> $sdir/TMP/chip/$line.chrom &
  done < "${EXP}"

fi
wait
#preprocessing input

if [[ -d $sdir/TMP/input ]]; then
  printf "creating intermediate input files\n"
  while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line
    chrom=$line
    grep "$chrom\>" $cdir/*.bedGraph >> $sdir/TMP/input/$line.chrom &
  done < "${EXP}"
fi
wait

printf "DONE!\n"
##launch Rscript for each chromosome

#humanchrlist=( $EXP )
printf "$humanchrlist"\n

if [[ ! -z $cdir ]]; then
  while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line

               outnamefileR=$line
               printf "Running Rscript..... $sPath/PeakZCL.R -f $sdir/TMP/chip/$line.chrom -c $sdir/TMP/input/$line.chrom -d "$odir" -n "$normalization" -j "$clustering" -w "$decimating" -l "$length" -s "$scales" -g "$levth" -z "$enrichment" -a "$area" -k "$foldchange" -m "$allpeaks" -o "$outnamefileR" \n" &
               Rscript $sPath/PeakZCL.R -f $sdir/TMP/chip/$line.chrom -c $sdir/TMP/input/$line.chrom -d "$odir" -n "$normalization" -w "$decimating" -l "$length" -s "$scales" -g "$levth" -z "$enrichment" -j "$clustering" -a "$area" -k "$foldchange" -m "$allpeaks" -o "$outnamefileR" &
               NPROC=$(($NPROC+1))
               if [ "$NPROC" -ge "$cores" ]; then
                 wait
                 NPROC=0
               fi

  done < "${EXP}"
else
    while IFS='' read -r line || [[ -n "$line" ]]; do
              outnamefileR=$line
              printf "Running Rscript ..... $sPath/PeakZCL.R -f $sdir/TMP/chip/$line.chrom -d "$odir" -n "$normalization" -j "$clustering" - "$decimating" -l "$length" -s "$scales" -g "$levth" -z "$enrichment" -a "$area" -k "$foldchange" -m "$allpeaks" -o "$outnamefileR" \n" &
              Rscript $sPath/PeakZCL.R -f $sdir/TMP/chip/$line.chrom -d "$odir" -n "$normalization" -w "$decimating" -j "$clustering" -l "$length" -s "$scales" -g "$levth" -z "$enrichment" -a "$area" -k "$foldchange" -m "$allpeaks" -o "$outnamefileR" &
              NPROC=$(($NPROC+1))
              if [ "$NPROC" -ge "$cores" ]; then
                wait
                NPROC=0
              fi
  done < "${EXP}"
fi

wait

printf "concatenate output files\n"
#concatenate files
cat $odir*.bed >> $odir$outnamefile.bed
sort -k1,1V -k2,2g -n $odir$outnamefile.bed >> $odir$outnamefile"_"$levth"_"$area"_"$clustering"_"$foldchange"_"$normalization.bed
#rm $odir$outnamefile.bed
#Rscript $sPath/process_Output_bed.R $odir$outnamefile.bed

peaks= wc -l $odir$outnamefile"_"$levth"_"$area"_"$clustering"_"$foldchange"_"$normalization.bed
#rm $odir$outnamefile"_ordered".bed
printf "$peaks peaks found\n"

if [[ ! -z $annotfile ]]; then
  Rscript $sPath/Annotation.R $annotfile $odir$outnamefile"_"$levth"_"$area"_"$clustering"_"$foldchange"_"$normalization.bed $odir$outnamefile"_"$levth"_"$area"_"$clustering"_"$foldchange"_"$normalization"_annot.txt"

fi




#rm -rf $sdir/TMP

day=`date +"%d/%m/%Y"`
hour=`date +"%H:%M"`
printf "finished at $day $hour !\n"
