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
| MODEL | Defines the choice of model considered during MC simulation. <br> *(Ising / XYZ)* | Ising |
| MCS | Monte Carlo steps per spins (MCS). <br> *[Total MCS (TMCS) steps for MC simulation, Total MCS steps for equilibration (TMCS_EQ) ]* | 5000 3000 |
| TEMP | Sets the temperature range. <br> *Note: Simulation will start from higher temperature.* <br> (T_start, T_end, T_interval) | 100 10 5 |
| SPIN | Spin values of each ion present in the `structure.vasp` file. <br> *Note: For non-magnetic ions, assign 0.0* |  |
| SPECIES | Symbol of species for which MC simulation is intended |  |
| SC | Size of supercell along x, y, z direction *(X, Y, Z)* | 2 2 2 |
| STG | Logic for staggered magnetization. If `.TRUE.` then `ETHER` requires the mandatory `staggered` file. <br> (Logic) | `.FALSE.` <br> cᵢ = 1 |
| BC | Boundary conditions (c: Closed, o: Open) along x, y, z direction *(X, Y, Z)* | c c c |
| REPEAT | Repeat MC simulations *m* times for better averaging results | 10 |
| ANGLE | Updates direction of spin vector. <br><br> If `.TRUE.`: <br> Sₓ = √(1 − S_z²) cosφ <br> Sᵧ = √(1 − S_z²) sinφ <br> −1 ≤ S_z ≤ +1 <br> φ = φ_old + φ'δ (0 ≤ δ ≤ 1) <br><br> Otherwise, spins follow George Marsaglia rule. | `.FALSE.` 5 |
| COA | Calculate observables at each defined MCS | 10 |
| ZEEMAN | Zeeman logic with external magnetic field <br> *(logic, Hₓ, Hᵧ, H_z)* | `.FALSE.` 0 0 0 |
| G_FACTOR | g factor value | 2 |
| SIA | Single ion anisotropy <br> (logic) | `.FALSE.` |
| PARA | Parameter logic with J_para (in meV) | `.FALSE.` 1.0 |
| OVRR | Over relaxation method applied after each ovrr_MCS steps <br> *(logic, over_para, ovrr_MCS)* | `.FALSE.` 0 0 |
| SEED | Seed for random number generation | 1992 |
| NBDFC | Neighbourhood finding criteria | 10⁻⁵ |
| PT |  Parallel Tempering <br> (logic, exchange_interval) | .FALSE. 11 |
| Checkerboard | Choice of STG lists in Checkerboard fashion. <br> (logic) | .FALSE. |
| Random_Algo | Choice of Random number generator algorithms. <br> (Random_Algo = mt19937/xoshiro256) | mt19937 |
| Scaled | Scaled the spins from classical (S) to qauntum level √S(S + 1) <br> _see:_ [arXiv:2602.16501](https://doi.org/10.48550/arXiv.2602.16501 "Click to Read") <br> (logic) | .False. |
| TB2J | Logic to allow ETHER to read TB2J (_see:_ [TB2J’s documentation!](https://tb2j.readthedocs.io/en/stable/index.html "Click to Read"))  output file 'exchange.out' (must be present in the working directory) to extract the interactios (J_iso, J_ani, DMI) for the ion pairs with bond length ≤ maxbond. <br> _Note: Currently Hamiltonian part containing J_iso, J_ani, and DMI are considered directly from TB2J documentation. ALong with these, one can also add Hamiltonian terms such as SIA and Zeemans terms by providing the details into the input files._ <br> (logic, maxbond) | .False. 0 |
| Home | Restricting the extractions (_from TB2J output file: exchange.out_) of J and other interaction parameters for the ion pairs present only in the Home cell (R: [0, 0, 0])  <br> (logic) | .True. |


## Monte Carlo with the TB2J Hamiltonian

Parameters $`J_{ij}`$, $`\mathbf{J}^{ani}_{ij}`$ and $`\vec{D}_{ij}`$ are taken directly from TB2J (`exchange.out`, in meV), without any conversion. Spins are unit vectors, $`\lvert \vec{S}_i \rvert = 1`$.

TB2J Hamiltonian (sum over all ordered pairs, so each bond appears twice, as $`(i,j)`$ and $`(j,i)`$):

```math
E = -\sum_i K_i \left(\vec{S}_i \cdot \vec{e}_i\right)^2 - \sum_{i,j} \left[ J_{ij} \vec{S}_i \cdot \vec{S}_j + \vec{S}_i \mathbf{J}^{ani}_{ij} \vec{S}_j + \vec{D}_{ij} \cdot \left( \vec{S}_i \times \vec{S}_j \right) \right]
```

### 1. Site energy of spin $`i`$ (for MC updates)

```math
E_i = -K_i \left(\vec{S}_i \cdot \vec{e}_i\right)^2 - \sum_{j} \left[ 2J_{ij} \vec{S}_i \cdot \vec{S}_j + 2 \vec{S}_i \mathbf{J}^{ani}_{ij} \vec{S}_j + 2 \vec{D}_{ij} \cdot \left( \vec{S}_i \times \vec{S}_j \right) \right]
```

Metropolis update $`\vec{S}_i \rightarrow \vec{S}_i'`$:

```math
\Delta E = E_i(\vec{S}_i') - E_i(\vec{S}_i), \qquad P_{acc} = \min\left(1, e^{-\Delta E / k_B T}\right)
```

### 2. Total energy

Sum over all ordered pairs (each bond appears twice, as $`(i,j)`$ and $`(j,i)`$):

```math
E_{tot} = -\sum_i K_i \left(\vec{S}_i \cdot \vec{e}_i\right)^2 - \sum_i \sum_{j} \left[ J_{ij} \vec{S}_i \cdot \vec{S}_j + \vec{S}_i \mathbf{J}^{ani}_{ij} \vec{S}_j + \vec{D}_{ij} \cdot \left( \vec{S}_i \times \vec{S}_j \right) \right]
```


