#!/bin/bash

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
   -d   database (GOBP, GOcc, etc....)
   -k   datos (TCGA, GTEX, CCLE...)
   -p   procesors
EOF
}

#Defaults --

while getopts "f:d:p:" OPTION
do
	case $OPTION in
	f) sdir=$OPTARG
	;;
	d) database=$OPTARG
	;;
    k) datos=$OPTARG
    ;;)
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

directory=$sdir"chunks_"$database"/"
cd $directory
Files=$(find . -maxdepth 1 -iname "*.Rdata")
printf '%s\n' "${Files[@]}" > $sdir"/chunk_file_list_$database.txt"
if [[ ! -z $sdir'OUTPUT'$datos"_"$database"/" ]]; then
    mkdir $sdir'OUTPUT'$datos"_"$database"/"
fi


#fi
export LANG=C
export LC_ALL=C

EXP=$sdir"chunk_file_list.txt"

while IFS='' read -r line || [[ -n "$line" ]]; do
              outnamefileR=$line
              printf "Running Rscript.....$sPath/uPEs_automatizer_step2.R  "$sdir" "$database" "$datos" "$outnamefileR"  \n" &
              Rscript $sPath/uPEs_automatizer_step2.R "$sdir" "$database" "$datos" "$outnamefileR" &
              NPROC=$(($NPROC+1))
              if [ "$NPROC" -ge "$cores" ]; then
                wait
                NPROC=0
              fi
  done < "${EXP}"
fi

wait

#rm -rf $sdir/TMP

day=`date +"%d/%m/%Y"`
hour=`date +"%H:%M"`
printf "finished at $day $hour !\n"
