#!/bin/bash
loc="$PWD/ether"
mpirun="nohup mpirun -np"
out="out"


export OMP_NUM_THREADS=2
number_of_processor=6
ulimit -v unlimited
ulimit -s unlimited


echo "$mpirun $number_of_processor \"$loc\" < /dev/null > $out &" > run.sh
chmod +x run.sh
echo 'Job is RUNNING'
./run.sh
sleep 1
rm -f run.sh

# For runnning the Ether, executable ./ether should be present in the woring directory

# Please note: Ether can also run via simple set of lines. Just remove the commented sign '#' 
# and paste into the other file [say 'run_job.sh'] and make it executable by typing chmod +x run_job.sh in the terminal and execute it by ./run_job.sh


#export OMP_NUM_THREADS=2
# ulimit -v unlimited
# ulimit -s unlimited
# mpirun -np 6 ./ether
