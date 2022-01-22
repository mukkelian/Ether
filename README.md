# Ether
Based on the classical lattice model (Heisenberg, XY, XYZ, etc.), code Ether has been developed to study the thermodynamics of ANY CRYSTAL SYSTEM by performing the basic Monte Carlo methods. Metropolis algorithm has been used to equate all the observables.

Code Ether solves the given Hamiltonian
H = J.SiSj - g*mb*MSi	-----(1)

or,
H = Jxx.(SxiSxj) + Jyy.(SyiSyj) + Jzz.(SziSzj) - g.mb.(MxSx + MySy + MzSz)	-----(2)

where, 

a is anisotrpy factor, a--> (a1, a2, a3) :: (on Jxx, on Jyy, on Jzz)

J is defined as exchang energy (unit meV) splited into Jxx, Jyy, Jzz components.
Jxx = a1.J_value, Jyy = a2.J_value, Jzz = a3.J_value. **J_value** is a input variable, can be supplied in '**j_exchange**' input file

S is spin vector and has Sx, Sy, Sz components.

M is magnetic field vector (unit Tesla) and has Mx, My, Mz components as well.

g is g-factor

mb is Bhor magneton (ev/T)

**Note::** Always supplied 'J_value' values in unit of meV (milli electron volts). For example for 0.001 eV provide this value in Ether code as 1, code will automatically convert it into the meV. For classical Heisenberg hamiltonian, kindly put all a1, a2, a3 equal to 1 which represent Jxx=Jyy=Jzz.

Compile this **ether.f90** by any FORTRAN compiler, for our case we have choosen **gfortran**

For the compilation, type the command same like given below into the terminal

 gfortran ether.f90 -o ether

executable 'ether' will be created for running this program.

# 1. About 'input' file

In input file there are basic informations which you need to provide before running the program which are given below,

I		---> Ising/Heisenberg (I/M) model

20000		---> MC steps

10000		---> equilibration steps

100 2 2		---> Temp(K) ==> final, initial, interval

3		---> no. of species present in the structure.vasp file 

0.5 0 0		---> Magnetic moment for species1, species2, species3.. and so on..

1		---> no. of species to include

Ce		---> species symbols (repeat in same line for many species)

4 4 4		---> Supercell size

.F.		---> for staggered  magn. (optional)

c c c		---> boundary CLOSED/OPEN (c/o)

2		---> for sampling (used to calculate the statistical error); NOTE :: you can increase it for correct result

50		---> least phi angle (used to achieve the equilibrium spin state, phi angle will vary within this range) 

50		---> MC step interval to calculate observables

.F. 5 5 0	---> Magnetic field (logic, Mx, My, Mz)

2		---> g_factor 

# 2. About 'j_exchange' file

4		---> no. of distinct nbd(N) 

Ce ce 1.0 1 1 1 4.16970	(for N=1)	! ith ion, jth ion, J_value, anisotropy(a1, a2, a3), bond length.

ce Ce 1.0 1 1 1 7.22213	(for N=2)

Ce Ce 1.0 1 1 1 4.08267	(for N=3)

ce ce 1.0 1 1 1 7.57543	(for N=4)

#note:	small/large caps are automatically considered into the program. For anisotropy see (a1, a2, a3) in equation no.2. a can be set from 0 <--> 1 and will considered as like 0 <--> 100 % 

# 3. About structure.vasp file

'structure.vasp' is an input structure file used to perform the Density Functional Theory (DFT) calculation using VASP program ( https://www.vasp.at/ ). For the making structure file use VESTA tool (https://jp-minerals.org/vesta/en/download.html) and export it into the file name 'structure.vasp'.
  # Please always export the file into cartesian coordinates.
  
Contact me if somebody face problem in understanding the input files.

Mukesh Kumar Sharma
email ID:: msharma1@ph.iitr.ac.in

If you find this code helpful in your research work kindly cite this code. New collaborations will be welcomed.

# 4. Output section
Along with all data files, there are two directories (_data, _spin) which will be generated during execution. 

# 4.1 About _spin directory
It contains the final spin configuration files for each temperature. File is named as "sp***K.dat", K is referred as Kelvin.

To visualize these files, please install xcrysden software (http://www.xcrysden.org/Download.html). Please follow the given step to visualize the spin configuration as well as lattice distribution.

Step 1.
open xcrysden software by typing **./xcrysden** in terminal.

Step 2.
go to the file --> Open structure --> Open XSF (Xcrysdgen structure file)
Step 3. Select the desired sp***K.dat file

Step 4.
At this point, you can see the lattice arrangement. Now, press 'f' button then the 'Shift + f' button combination. After pressing the 'Shift + f' you will have another dialogue box, in that there is an option to change the 'Length factor'. Change it by 4 or 5 or any number you want. You will have the spin configuration for the selected temperature.

# 4.2 About _data directory

In this directory data files related to energy, magnetization, nbd (neighborhood), and initial spin configuration are present.

For the basic plot where information of specific heat, susceptibility, magnetization, energy/site, just type 'gnuplot graph.sh' in the terminal, the output file named as fig1.png will be generated.

please always check the nbd.dat file to make sure that code is taking the right nbd's for the selected lattice point. To do this open the initial_spin_conf.dat file into the Xcrysden and this 'nbd.dat' file in gedit/Notepad editor and check and compare the nbd information present in the nbd.dat file. 
(by clicking the Atom info. in Xcrysden you will have the lattice point ID on the screen, by this you can check the ID's of nbd and compare it from the 'nbd.dat' file)
