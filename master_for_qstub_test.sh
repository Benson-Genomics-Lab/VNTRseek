#!/bin/bash

if [ -z "$ARG1" ];
then
     echo "Error: dbsufix is not passed!"
     exit 1
fi

if [ -z "$ARG2" ];
then
     echo "Error: nprocesses is not passed!"
     exit 1
fi

if [ -z "$ARG3" ];
then
     echo "Error: path is not passed!"
     exit 1
fi

perl ${ARG3}/vntrseek.pl 0 19 --dbsuffix $ARG1 --nprocesses $ARG2 &
wait
