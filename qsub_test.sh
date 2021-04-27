#!/bin/bash
if (($# != 2));
then
     echo "Error: Please pass DBSUFIX and NPROCESSORS to this script. Aborting!"
     exit 1
fi
WORKD="$( dirname "${BASH_SOURCE[0]}" )"
qsub -cwd -l mem_total=64G -l h_rt=40:00:00 -pe omp $2 -v ARG1=$1,ARG2=$2,ARG3=$WORKD ${WORKD}/master_for_qstub_test.sh
