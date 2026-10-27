#!/bin/bash
loc="$PWD/ether"
mpirun="setsid nohup mpirun -np"
out="out"

number_of_processor=8i		# MPI (mpirun)
export OMP_NUM_THREADS=2	# OpenMP (threads)
export OMP_STACKSIZE=512M

mpi_opts="--use-hwthread-cpus --map-by slot:PE=$OMP_NUM_THREADS --bind-to hwthread"
ulimit -v unlimited
ulimit -s unlimited


echo "$mpirun $number_of_processor $mpi_opts \"$loc\" < /dev/null > $out 2>&1 &" > run.sh
echo 'echo $! > job.pid' >> run.sh
chmod +x run.sh
./run.sh
sleep 0
rm -f run.sh
echo "Job is RUNNING (PID $(cat job.pid))"
echo "To stop the whole job:  kill -- -\$(cat job.pid)"

