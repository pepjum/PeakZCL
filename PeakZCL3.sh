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
   -g   Genome reference file (required)
   -n   Normalization method ("SES","SDS", "N", Default: SES)
   -w   Decimating signal factor (Default: 50)
   -l   Length signal wavelet (Default: 200000)
   -s   Scales of the wavelet (Default: 40)
   -t   Zero-crossing lines leverage threshold (Default: 6)
   -z   Nature of enrichment. Options are ("histone" or "TF". Default: "TF")
   -b   Bp range between peaks to clustering purposes (Default: 1)
   -k   minimum width of the peak (Default:300)
   -x   maximum width of the peak (Default: 5000)
   -o   Output BED file name (Default: ZCL_peaks.final.bed)
   -v   Path to annotation database in .gtf format (Default: NULL)
   -m   Returns all peaks detected or not (VALUES: TRUE, FALSE. Default=FALSE)
   -p   Number of processors used by R scripts (Default: 1)
EOF
}

#Defaults --

normalization="SES"
decimating="50"
length="200000"
scales="40"
levth="10"
clustering="1"
klustering="1"
area="0"
foldchange="0.1"
cores="6"
sdir=""
cdir=""
genome=""
odir="NULL"
outnamefile="PeakZCL_peaks"
annotfile=""
enrichment="TF"
minpeak="300"
maxpeak="5000"
allpeaks="FALSE"

while getopts "f:c:n:w:l:s:g:t:z:j:b:a:k:x:o:v:y:m:p:" OPTION
do
	case $OPTION in
	f) sdir=$OPTARG
	;;
	c) cdir=$OPTARG
	;;
	n) normalization=$OPTARG
	;;
    g) genome=$OPTARG
    ;;
	w) decimating=$OPTARG
	;;
	l) length=$OPTARG
	;;
	s) scales=$OPTARG
	;;
	t) levth=$OPTARG
	;;
	z) enrichment=$OPTARG
	;;
	j) clustering=$OPTARG
    ;;
    b) klustering=$OPTARG
	;;
	a) area=$OPTARG
	;;
	k) minpeak=$OPTARG
	;;
    x) maxpeak=$OPTARG
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
if [[ -z $genome ]]
then
     usage
     exit 1
fi


# echo -e $sdir\n
# echo -e $cdir\n
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
# #wdir=$sdir'TMP/chip/'
# #kdir= $sdir'TMP/input/'
# #
if [[ $odir == "NULL" ]]; then
    mkdir $sdir'OUTPUT6/'
    outdir=$sdir'OUTPUT6/'
    find $outdir -maxdepth 1 -type f -delete

fi
export LANG=C
export LC_ALL=C

#preprocessing chip

if [[ ! -z $sdir ]]; then
  printf "creating coverage file chip file\n"
  genomeCoverageBed -d -split -ibam $sdir/*.bam -g $genome > $sdir/coverage_chip.bedGraph &


fi

if [[ ! -z $cdir ]]; then
  printf "creating coverage input file\n"
  genomeCoverageBed -d -split -ibam $cdir/*.bam -g $genome > $cdir/coverage_input.bedGraph &
fi

wait


if [[ ! -z $sdir ]]; then
  printf "creating intermediate chip files\n"

  while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line
    chrom=$line
    printf "$line\n"
     echo -e  "grep "$chrom\>" $sdir/*.bedGraph >> $sdir/TMP/chip/$line.chrom &"
    grep "$chrom\>" $sdir/*.bedGraph > $sdir/TMP/chip/$line.chrom &
  done < "${EXP}"

fi
wait
#preprocessing input

if [[ -d $sdir/TMP/input ]]; then
  printf "creating intermediate input files\n"
  while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line
    chrom=$line
    grep "$chrom\>" $cdir/*.bedGraph > $sdir/TMP/input/$line.chrom &
  done < "${EXP}"
fi
wait

printf "DONE!\n"
##launch Rscript for each chromosome

 if [[ ! -z $cdir ]]; then
   while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line

                outnamefileR=$line
                printf "Running Rscript..... $sPath/PeakZCL2.R -f $sdir/TMP/chip/$line.chrom -c $sdir/TMP/input/$line.chrom -n "$normalization" -j "$clustering" -w "$decimating" -l "$length" -s "$scales" -g "$levth" -z "$enrichment" -a "$area" -k "$minpeak" -x "$maxpeak" -m "$allpeaks" -o "$outnamefileR"  \n" &
                Rscript $sPath/PeakZCL2.R -f $sdir/TMP/chip/$line.chrom -c $sdir/TMP/input/$line.chrom -n "$normalization" -w "$decimating" -l "$length" -s "$scales" -g "$levth" -z "$enrichment" -j "$clustering" -a "$area" -k "$minpeak" -x "$maxpeak" -m "$allpeaks" -o "$outnamefileR"  &
                NPROC=$(($NPROC+1))
                if [ "$NPROC" -ge "$cores" ]; then
                  wait
                  NPROC=0
                fi

   done < "${EXP}"
 else
     while IFS='' read -r line || [[ -n "$line" ]]; do
               outnamefileR=$line
               printf "Running Rscript ..... $sPath/PeakZCL2.R -f $sdir/TMP/chip/$line.chrom -n "$normalization" -j "$clustering" -w "$decimating" -l "$length" -s "$scales" -g "$levth" -z "$enrichment" -a "$area" -k "$minpeak" -x "$maxpeak" -m "$allpeaks" -o "$outnamefileR"  \n" &
               Rscript $sPath/PeakZCL2.R -f $sdir/TMP/chip/$line.chrom -n "$normalization" -w "$decimating" -j "$clustering" -l "$length" -s "$scales" -g "$levth" -z "$enrichment" -a "$area" -k "$minpeak" -x "$maxpeak" -m "$allpeaks" -o "$outnamefileR"   &
               NPROC=$(($NPROC+1))
               if [ "$NPROC" -ge "$cores" ]; then
                 wait
                 NPROC=0
               fi
   done < "${EXP}"
 fi

 wait
 ###### launch bedtools multicov for each bed file

printf "Running bedtools multicov\n"
 if [[ ! -z $cdir ]]; then
   while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line

                outnamefileR=$line
                bedtools multicov -bams $sdir/*.bam -bed $sdir/TMP/chip/$line.bed > $sdir/TMP/chip/$line.txt &
                bedtools multicov -bams $cdir/*.bam -bed $sdir/TMP/chip/$line.bed > $sdir/TMP/input/$line.txt &
                NPROC=$(($NPROC+1))
                if [ "$NPROC" -ge "$cores" ]; then
                  wait
                  NPROC=0
                fi

   done < "${EXP}"
 else
     printf "Running bedtools multicov\n"
     while IFS='' read -r line || [[ -n "$line" ]]; do
               outnamefileR=$line
               bedtools multicov -bams $sdir/*.bam -bed $sdir/TMP/chip/$line.bed > $sdir/TMP/chip/$line.txt &

               NPROC=$(($NPROC+1))
               if [ "$NPROC" -ge "$cores" ]; then
                 wait
                 NPROC=0
               fi
   done < "${EXP}"
 fi

 wait

#### launch part 3

if [[ ! -z $cdir ]]; then
  while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line

               outnamefileR=$line

               printf "Running Rscript..... $sPath/Normalization.R -f $sdir/TMP/chip/$line.bed -c $sdir/TMP/chip/$line.txt -d $sdir/TMP/input/$line.txt -n $sdir/TMP/chip/$line.fc -p "$outdir" -o "$outnamefileR" -m "$allpeaks" -b "$klustering" -w "$sdir/TMP/chip/$line.chrom" \n" &
               Rscript $sPath/Normalization.R -f $sdir/TMP/chip/$line.bed -c $sdir/TMP/chip/$line.txt -d $sdir/TMP/input/$line.txt -n $sdir/TMP/chip/$line.fc -p "$outdir" -m "$allpeaks" -o "$outnamefileR" -b $klustering -w $sdir/TMP/chip/$line.chrom &

               NPROC=$(($NPROC+1))
               if [ "$NPROC" -ge "$cores" ]; then
                 wait
                 NPROC=0
               fi

  done < "${EXP}"
else
    while IFS='' read -r line || [[ -n "$line" ]]; do
              outnamefileR=$line

              printf "Running Rscript..... $sPath/Normalization.R -f $sdir/TMP/chip/$line.bed -c $sdir/TMP/chip/$line.txt -p "$outdir" -o "$outnamefileR" -b "$klustering" -w "$sdir/TMP/chip/$line.chrom" \n" &
              Rscript $sPath/Normalization.R -f $sdir/TMP/chip/$line.bed -c $sdir/TMP/chip/$line.txt -p "$outdir" -o "$outnamefileR" -b $klustering -w $sdir/TMP/chip/$line.chrom &


              NPROC=$(($NPROC+1))
              if [ "$NPROC" -ge "$cores" ]; then
                wait
                NPROC=0
              fi
  done < "${EXP}"
fi
wait
# #
# #
#
#
printf "concatenate output files\n"
# #concatenate files
cat $outdir*.Norm > $outdir$outnamefile.Norm
# # rm $outdir*.Norm
# # rm $outdir$outnamefile.bed
wait

printf "Running model.R ....\n"

Rscript $sPath/model.R $outdir$outnamefile.Norm
wait

printf "Running Poisson.R ....\n"
if [[ ! -z $cdir ]]; then

    control=TRUE
    Rscript $sPath/Poisson.R $outdir$outnamefile.model $allpeaks $klustering $control $maxpeak

    else
    control=FALSE
    Rscript $sPath/Poisson.R $outdir$outnamefile.model $allpeaks $klustering $control $maxpeak
fi
wait
#
sort -k1,1V -k2,2g $outdir$outnamefile"_output_poissonFDR.bed" > $outdir$outnamefile"_zCros_"$levth"_clus_"$clustering"_norm_"$normalization"_min_"$minpeak"_max_"$maxpeak.sorted_final.bed
peaks= wc -l $outdir$outnamefile"_zCros_"$levth"_clus_"$clustering"_norm_"$normalization"_min_"$minpeak"_max_"$maxpeak.sorted_final.bed
# # #rm $odir$outnamefile"_ordered".bed
# # printf "$peaks peaks found\n"
# #
# if [[ ! -z $annotfile ]]; then
#   Rscript $sPath/Annotation.R $outdir$outnamefile"_zCros_"$levth"_clus_"$clustering"_norm_"$normalization"_min_"$minpeak"_max_"$maxpeak.sorted_final.bed $outdir$outnamefile"_zCros_"$levth"_clus_"$clustering"_norm_"$normalization"_min_"$minpeak"_max_"$maxpeak.sorted_final_annot.txt
# fi
#



#rm -rf $sdir/TMP

day=`date +"%d/%m/%Y"`
hour=`date +"%H:%M"`
printf "finished at $day $hour !\n"
