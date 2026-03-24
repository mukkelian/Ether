# Efficient Tool for THermodynamics Exploration via Relaxations (ETHER)


Monte Carlo (MC) simulations are powerful computational tools for investigating thermodynamic behavior and validating analytical approaches in complex physical systems. Here we present ETHER (Efficient Tool for THermodynamics Exploration via Relaxations), an open-source MC simulation package, developed for studying temperature-dependent magnetic properties in spin systems. ETHER enables large-scale MC simulations of various spin systems to analyze phase transitions, critical behavior, and complex magnetic structures. The package constructs spin-lattice networks from standard structural input files and supports exchange interaction definitions through user-specified neighbor lists. It also provides tools for easy visualization and post-processing of simulation outputs. Our code has been benchmarked thoroughly against results reported in the literature for common representative magnetic systems. This user-friendly code offers researchers a versatile platform for exploring thermodynamic properties of complex magnetic systems.

*Contacts*
msharma1@ph.iitr.ac.in
mukelian92@gmail.com 


## Citation

If this repository helps your research or work, please consider citing our code.  
The citation information is available in [citation.bib](./citation.bib). 

## Overwiew of Input File Details ##

`$in.ether$` file contains the necessary input information based on which `ETHER` will build up the internal framework and initiate the MC simulation accordingly. Users can provide the following *Tags* for initiating the MC simulation:

| Tags     | Description | Default |
|----------|------------|---------|
| MODEL | Defines the choice of model considered during MC simulation. *(Ising / XYZ)* | Ising |
| MCS | Monte Carlo steps per spins (MCS). Total MCS for simulation and equilibration process | 5000 3000 |
| TEMP | Sets the temperature range. <br> *Note: Simulation will start from higher temperature.* <br> (T_start, T_end, T_interval) | 100 10 5 |
| SPIN | Spin values of each ion present in the `structure.vasp` file. <br> *Note: For non-magnetic ions, assign 0.0* |  |
| SPECIES | Symbol of species for which MC simulation is intended |  |
| SC | Size of supercell along x, y, z direction *(X, Y, Z)* | 2 2 2 |
| STG | Logic for staggered magnetization. If `.TRUE.` then `ETHER` requires the mandatory `staggered` file. <br> (Logic) | `.FALSE.` <br> cᵢ = 1 |
| BC | Boundary conditions (c: Closed, o: Open) along x, y, z direction *(X, Y, Z)* | c c c |
| REPEAT | Repeat MC simulations *m* times for better averaging results | 10 |
| ANGLE | Updates direction of spin vector. <br><br> If `.TRUE.`: <br> Sₓ = √(1 − S_z²) cosφ <br> Sᵧ = √(1 − S_z²) sinφ <br> −1 ≤ S_z ≤ +1 <br> φ = φ_old + φ'δ (0 ≤ δ ≤ 1) <br><br> Otherwise, spins follow George Marsaglia rule. | `.FALSE.` 5 |
| COA | Calculate observables at each defined MCS | 10 |
| ZEEMAN | Zeeman logic with external magnetic field *(logic, Hₓ, Hᵧ, H_z)* | `.FALSE.` 0 0 0 |
| G_FACTOR | g factor value | 2 |
| SIA | Single ion anisotropy logic | `.FALSE.` |
| PARA | Parameter logic with J_para (in meV) | `.FALSE.` 1.0 |
| OVRR | Over relaxation method applied after each ovrr_MCS steps *(logic, over_para, ovrr_MCS)* | `.FALSE.` 0 0 |
| SEED | Seed for random number generation | 1992 |
| NBDFC | Neighbourhood finding criteria | 10⁻⁵ |
