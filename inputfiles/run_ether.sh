#!/bin/bash
loc=$PWD/ether
loc1="nohup"
loc2="< /dev/null > out &"
echo "$loc1" "$loc" "$loc2" >> run.bat
chmod +x run.bat
echo 'Job is RUNNING'
OMP_NUM_THREADS=64
export OMP_NUM_THREADS
ulimit -s unlimited
./run.bat
rm *.bat
