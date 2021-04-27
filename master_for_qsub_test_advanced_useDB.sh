#!/bin/bash -l

if [ -z "$GID" ];
then
     echo "Error: GID is not passed!"
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

batchname=${WORKD}/../qsub_test_advanced_useDB.sh
# scriptname=${WORKD}/../master_for_qsub_test_advanced_useDB.sh
perlfile=${WORKD}/vntrseek.pl
MYSQLLOGIN="$(grep LOGIN "${WORKD}"/vs.cnf | cut -f2 -d'=')"
MYSQLPASS="$(grep PASS "${WORKD}"/vs.cnf | cut -f2 -d'=')"
MYSQLHOST="$(grep HOST "${WORKD}"/vs.cnf | cut -f2 -d'=')"

echo "$GID"

# ask for dbsuffix of run
runinfo=( $( (mysql -N -h "${MYSQLHOST}" -u "${MYSQLLOGIN}" -p"${MYSQLPASS}" vntrseektrack) <<SQL
LOCK TABLES ${dbtable} READ;
SELECT name, dbsuffix FROM ${dbtable} WHERE gid = $GID;
SQL
) )

DBSUFFIX=${runinfo[1]}
echo "$DBSUFFIX"

STEPS=""

handle_error() {
  runnext=$(perl "$perlfile" 99 --dbsuffix "$DBSUFFIX")
  laststep=$(($runnext - 1))
  echo "Oh noes, something went wrong (got SIGERR)! Setting run to status 2 (error)"
  mysql -N -h "${MYSQLHOST}" -u "${MYSQLLOGIN}" -p"${MYSQLPASS}" vntrseektrack <<SQL
LOCK TABLES ${dbtable} WRITE;
UPDATE ${dbtable} SET status=2, isrunning=0, laststep=$laststep WHERE gid=$GID;
SQL
  echo "Job encountered error at step $laststep"
  exit 1
}

handle_kill() {
  # A kill signal was sent.
  runnext=$(perl "$perlfile" 99 --dbsuffix "$DBSUFFIX")
  laststep=$(($runnext - 1))
  echo "Oh noes, something went wrong (got SIGKILL)! Setting run to status 2 (error)"
  mysql -N -h "${MYSQLHOST}" -u "${MYSQLLOGIN}" -p"${MYSQLPASS}" vntrseektrack <<SQL
LOCK TABLES ${dbtable} WRITE;
UPDATE ${dbtable} SET status=2, isrunning=0, laststep=$laststep WHERE gid=$GID;
SQL
  echo "Job encountered error at step $laststep"
  exit 1
}

handle_success() {
  STATUS=""
  runnext=$(perl "$perlfile" 99 --dbsuffix "$DBSUFFIX")
  laststep=$(($runnext - 1))
  if [ $laststep -eq 19 ]; then
    echo "Last step! Setting run to status 1 (complete)"
    STATUS="status=1,"
  fi
  mysql -N -h "${MYSQLHOST}" -u "${MYSQLLOGIN}" -p"${MYSQLPASS}" vntrseektrack <<SQL
LOCK TABLES ${dbtable} WRITE;
UPDATE ${dbtable} SET $STATUS isrunning=0, laststep=$laststep WHERE gid=$GID;
SQL
  bash "$batchname" "$NPROCS" "${dbtable}"
}

trap handle_error ERR
trap handle_kill KILL

if [ "$runnext" -gt 19 ]; then
  echo "Nothing to do, done!"
  exit 0
elif [ "$runnext" -gt 17 ]; then
  echo "Doing steps 18 and 19!"       # single
  STEPS="18 19"
elif [ "$runnext" -gt 15 ]; then
  echo "Doing steps 16 and 17!"       # single
  STEPS="16 17"
elif [ "$runnext" -gt 14 ]; then
  echo "Doing step 15!"               # multiple
  STEPS="15"
elif [ "$runnext" -gt 13 ]; then
  echo "Doing step 14!"               # single
  STEPS="14"
elif [ "$runnext" -gt 12 ]; then
  echo "Doing step 13!"               # multiple
  STEPS="13"
elif [ "$runnext" -gt 10 ]; then
  echo "Doing steps 11-12!"           # single
  STEPS="11 12"
elif [ "$runnext" -gt 9 ]; then
  echo "Doing step 10!"               # multiple
  STEPS="10"
elif [ "$runnext" -gt 4 ]; then
  echo "Doing steps 5-9!"             # single
  STEPS="5 9"
elif [ "$runnext" -gt 3 ]; then
  echo "Doing step 4!"                # multiple
  STEPS="4"
elif [ "$runnext" -gt 1 ]; then
  echo "Doing steps 2-3!"             # single
  STEPS="2 3"
else
  echo "Doing steps 0-1!"             # multiple
  STEPS="0 1"
fi

if [ "$STEPS" ]; then
  # Can't use "" around $STEPS; not seen as two space-separated args
  perl "${perlfile}" ${STEPS} --dbsuffix "${DBSUFFIX}" --nprocesses ${NPROCS} &
  wait "$!"
  if [ $? -ne 0 ]; then
    handle_error "$DBSUFFIX"
  fi
  handle_success "$DBSUFFIX"
else
  echo "STEPS is not defined!"
  handle_error "$DBSUFFIX"
fi
