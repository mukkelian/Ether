loc="$PWD/ether"
mpirun="nohup mpirun -np"
number_of_processor=6
out="out"
echo "$mpirun $number_of_processor \"$loc\" < /dev/null > $out &" > run.sh
chmod +x run.sh
echo 'Job is RUNNING'
export OMP_NUM_THREADS=2
ulimit -s unlimited
./run.sh
sleep 1
rm -f run.sh
