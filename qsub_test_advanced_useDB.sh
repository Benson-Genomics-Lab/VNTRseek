#!/bin/bash -l

if (($# != 2));
then
     echo "Usage: $0 <NPROCS> <DBTABLE>. Aborting!"
     exit 1
fi

nprocs=$1
dbtable="$2"

handle_error() {
  if [ "$GID" ]; then
    mysql -N -h "${MYSQLHOST}" -u "${MYSQLLOGIN}" -p"${MYSQLPASS}" vntrseektrack <<SQL
LOCK TABLES ${dbtable} AS setRunning WRITE;
UPDATE ${dbtable} AS setRunning SET isrunning = 0 WHERE gid = $GID;
UNLOCK TABLES;
SQL
  fi
}

trap handle_error ERR

SCRIPTD="$( dirname "${BASH_SOURCE[0]}" )"
WORKD="${SCRIPTD}"
scriptname="${SCRIPTD}/master_for_qsub_test_advanced_useDB.sh"
# perlfile=${WORKD}/vntrseek.pl
MYSQLLOGIN="$(grep LOGIN "${WORKD}"/vs.cnf | cut -f2 -d'=')"
MYSQLPASS="$(grep PASS "${WORKD}"/vs.cnf | cut -f2 -d'=')"
MYSQLHOST="$(grep HOST "${WORKD}"/vs.cnf | cut -f2 -d'=')"

# perl $perlfile 99 --dbsuffix $1 --nprocesses $2   # ask for what step needs to be run next
# Lock table, ask for what genome and step needs to be run next,
# set it to running, and unlock table.
runinfo=( $( (mysql -N -h "${MYSQLHOST}" -u "${MYSQLLOGIN}" -p"${MYSQLPASS}" vntrseektrack) <<SQL
LOCK TABLES ${dbtable} AS setRunning WRITE, ${dbtable} AS getGid1 READ, ${dbtable} AS getGid2 READ;
SELECT gid,laststep FROM ${dbtable} AS getGid1 WHERE status = 0 AND isrunning = 0 LIMIT 1;
UPDATE ${dbtable} AS setRunning SET isrunning = 1
  WHERE gid IN (SELECT gid
    FROM (SELECT gid FROM ${dbtable} AS getGid2 
      WHERE status = 0 AND isrunning = 0
      LIMIT 1)
    AS t);
UNLOCK TABLES;
SQL
) )

if [ -z "${runinfo[0]}" ]; then
  echo "No results from MySQL (no more samples to run?). Stopping..."
  exit 0
fi

runnext=$(( ${runinfo[1]}+1 )) # add 1 to the last step run. Bash arrays are 0-indexed.

echo "${runinfo[@]}"
GID=${runinfo[0]}

case  $runnext  in
  0|4|10|13|15) 				# multiple
    echo "runnext state ($runnext), running in MULTI mode!"
    #qsub -cwd -m a -l mem_total=64G -l h_rt=50:00:00 -pe omp "$nprocs" -v GID=$GID,NPROCS="$nprocs",WORKD="$WORKD",runnext=$runnext,dbtable=$dbtable "$scriptname"
    qsub -cwd -m a -l mem_total=64G -l h_rt=40:00:00 -pe omp "$nprocs" -v GID=$GID,NPROCS="$nprocs",WORKD="$WORKD",runnext=$runnext,dbtable=$dbtable "$scriptname"
    #qsub -q lbi,bioinfo -cwd -m a -l mem_total=64G -l h_rt=50:00:00 -pe omp "$nprocs" -v GID=$GID,NPROCS="$nprocs",WORKD="$WORKD",runnext=$runnext,dbtable=$dbtable "$scriptname"
    ;;
  2|3|5|6|7|8|9|11|12|14|16|17|18|19)      	# single
    echo "runnext state ($runnext), running in SINGLE mode!"
    #qsub -cwd -m a -l h_rt=30:00:00 -v GID=$GID,NPROCS="$nprocs",WORKD="$WORKD",runnext=$runnext,dbtable=$dbtable "$scriptname"
    qsub -cwd -m a -l h_rt=30:00:00 -v GID=$GID,NPROCS="$nprocs",WORKD="$WORKD",runnext=$runnext,dbtable=$dbtable "$scriptname"
    #qsub -q lbi,bioinfo -cwd -m a -l h_rt=20:00:00 -v GID=$GID,NPROCS="$nprocs",WORKD="$WORKD",runnext=$runnext,dbtable=$dbtable "$scriptname"
    ;;
  1)
    echo "runnext state ($runnext), running in MULTI mode!"
    qsub -cwd -m a -l mem_total=64G -l h_rt=80:00:00 -pe omp "$nprocs" -v GID=$GID,NPROCS="$nprocs",WORKD="$WORKD",runnext=$runnext,dbtable=$dbtable "$scriptname"
    ;;
  20)
    echo "No more qsub jobs needed ($runnext), done!"
    ;;
  *)
    echo "Illegal runnext state ($runnext)!"
    ;;
esac 

