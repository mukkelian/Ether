export OMP_NUM_THREADS=2
ulimit -v unlimited
ulimit -s unlimited
mpirun -np 6 ./ether
