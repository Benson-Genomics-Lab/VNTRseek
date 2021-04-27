#!/bin/bash -l

if [ -z "$DBSUFFIX" ];
then
     echo "Error: dbsuffix is not passed!"
     exit 1
fi

if [ -z "$NPROCS" ];
then
     echo "Error: nprocesses is not passed!"
     exit 1
fi

if [ -z "$WORKD" ];
then
     echo "Error: path is not passed!"
     exit 1
fi

if [ -z "$runnext" ];
then
     echo "Error: next step is not passed!"
     exit 1
fi

batchname=${WORKD}/qsub_test_advanced.sh
perlfile=${WORKD}/vntrseek.pl

runnext=$(perl "$perlfile" 99 --dbsuffix "$DBSUFFIX") # ask for what step needs to be run next

STEPS=()

handle_error() {
  runnext=$(perl "$perlfile" 99 --dbsuffix "$DBSUFFIX")
  laststep=$(("$runnext" - 1))
  echo "Oh noes, something went wrong (got SIGERR)! Setting run to status 2 (error)"
  echo "Job encountered error at step $laststep"
  exit 1
}

handle_success() {
  # STATUS=""
  runnext=$(perl "$perlfile" 99 --dbsuffix "$DBSUFFIX")
  laststep=$((runnext - 1))
  if [ $laststep -eq 19 ]; then
    echo "Last step!"
    exit
  fi
  bash "$batchname" "$NPROCS" "${DBSUFFIX}"
}

trap handle_error ERR

if [ "$runnext" -gt 19 ]; then
  echo "Nothing to do, done!"
  exit 0
elif [ "$runnext" -gt 17 ]; then
  echo "Doing steps 18 and 19!"       # single
  STEPS=(18 19)
elif [ "$runnext" -gt 15 ]; then
  echo "Doing steps 16 and 17!"       # single
  STEPS=(16 17)
elif [ "$runnext" -gt 14 ]; then
  echo "Doing step 15!"               # multiple
  STEPS=(15)
elif [ "$runnext" -gt 13 ]; then
  echo "Doing step 14!"               # single
  STEPS=(14)
elif [ "$runnext" -gt 12 ]; then
  echo "Doing step 13!"               # multiple
  STEPS=(13)
elif [ "$runnext" -gt 10 ]; then
  echo "Doing steps 11-12!"           # single
  STEPS=(11 12)
elif [ "$runnext" -gt 9 ]; then
  echo "Doing step 10!"               # multiple
  STEPS=(10)
elif [ "$runnext" -gt 4 ]; then
  echo "Doing steps 5-9!"             # single
  STEPS=(5 9)
elif [ "$runnext" -gt 3 ]; then
  echo "Doing step 4!"                # multiple
  STEPS=(4)
elif [ "$runnext" -gt 1 ]; then
  echo "Doing steps 2-3!"             # single
  STEPS=(2 3)
else
  echo "Doing steps 0-1!"             # multiple
  STEPS=(0 1)
fi

if [ "${STEPS[@]}" ]; then
  perl "${perlfile}" "${STEPS[@]}" --dbsuffix "$DBSUFFIX" --nprocesses "$NPROCS" &
  wait "$!"
  if ! wait "$!"; then
    handle_error "$DBSUFFIX"
  fi
  handle_success "$DBSUFFIX"
else
  echo "STEPS is not defined!"
  handle_error "$DBSUFFIX"
fi