# Ether
Based on the classical lattice model (Heisenberg, XYZ, etc.), code Ether has been developed to study the thermodynamics of ANY CRYSTAL SYSTEM by performing the basic Monte Carlo methods. Metropolis algorithm has been used to equate all the observables

# 1. About 'input' file

In input file there are basic information which you need to provide before running the program as given below,

20000		  ! MC steps
10000		  ! equilibration steps
100 2 2		! Temp(K) ==> final, initial, interval
3		      ! no. of species present in structure.vasp file 
0.5 0 0   ! Magnetic moment for species1,species2, species3.. and so on..
1		      ! no. of species to include
Ce		    ! species symbols (repeat in same line for many species)
4 4 4		  ! Supercell size
.F.		    ! for staggered  magn. (optional)
c c c		  ! boundary CLOSED/OPEN (c/o)
2		      ! for sampling (used to calculate the statistical error); NOTE :: you can increase it for correct result
50		    ! least phi angle (used to achieve the equilibrium spin state, phi angle will vary within this range) 
50		    ! MC step interval to calculate observables
.F. 5 5 0	! Magnetic field (logic, Mx, My, Mz)
2		      ! g_factor 

# 2. About 'j_exchange_input' file

1			          ! no. of distinct nbd
4.16970         ! bond length for ith distinct nbd (d_nbd1, d_nbd2 ...)
1			          ! no. of lines present below
1 1 1 1.0 1 1 1	! (ith nbd, j-th species, k-th species, J_(k, j), Jxx, Jyy, Jzz); note:: provide Jxx/Jyy/Jzz (0-1) consider as 0% - 100%

# 3. About structure.vasp file
'structure.vasp' is a structure file used in performing the Density Functional Theory (DFT) by VASP ( https://www.vasp.at/ ). For the making structure file use VESTA tool (https://jp-minerals.org/vesta/en/download.html) and export it into the file name 'structure.vasp'.
  # Don't forget to export the file into cartesian coordinates not in fractinal coordinates
  
Contact me if somebody face problem in understanding the input files.
Mukesh Kumar Sharma
email ID:: msharma1@ph.iitr.ac.in
If you find this code helpful in your research work please cite this code.
