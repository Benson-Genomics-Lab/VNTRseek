#!/bin/bash -l

if (($# != 2));
then
     echo "Error: Please pass NPROCESSORS and DBSUFFIX to this script. Aborting!"
     exit 1
fi

WORKD="$( dirname "${BASH_SOURCE[0]}" )"
scriptname="$WORKD/master_for_qstub_test_advanced.sh"
perlfile="$WORKD/vntrseek.pl"

nprocs=$1
dbsuffix="$2"

runnext=$(perl "$perlfile" 99 --dbsuffix "$dbsuffix")   # ask for what step needs to be run next

case  $runnext  in
  0|1|4|10|13|15) 				# multiple
    echo "runnext state ($runnext), running in MULTI mode!"
    qsub -cwd -N vntrseek_"$dbsuffix"_step"$runnext" -l mem_total=64G -l h_rt=20:00:00 -pe omp "$nprocs" -v DBSUFFIX="$dbsuffix",NPROCS="$nprocs",WORKD="$WORKD",runnext="$runnext" "$scriptname"
    ;;
  2|3|5|6|7|8|9|11|12|14|16|17|18|19)      	# single
    echo "runnext state ($runnext), running in SINGLE mode!"
    qsub -cwd -N "vntrseek_${dbsuffix}_step$runnext" -l h_rt=20:00:00 -v DBSUFFIX="$dbsuffix",NPROCS="$nprocs",WORKD="$WORKD",runnext="$runnext" "$scriptname"
    ;;
  20)
    echo "No more qsub jobs needed ($runnext), done!"
    ;;
  *)
    echo "Illegal runnext state ($runnext)!"
    ;;
esac 

