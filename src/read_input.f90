! Ether, a Monte Carlo simulation program, impowers the users to study 
! the thermodynamical properties of spins arranged in any complex 
! lattice network.

! Copyright (C) 2021  Mukesh Kumar Sharma (msharma1@ph.iitr.ac.in)

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program; if not, see https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html

        subroutine read_input

	use init

	implicit none

	integer :: i
	logical :: file_present
		
	inquire(file='in.ether', exist=file_present)
	if(.not.file_present) then
	if (rank == 0) then
		write(6, *) "==> Required file 'in.ether' is not present"
		write(6, *) "	 STOPPING now"
		write(6, *) ""
	end if
		stop
	end if

        open(unit=0, file='in.ether', status='old', action='read')

	read(0, *) model				! Ising/XYZ model

	! CONVERTING STRING VARIBLES INTO LOWER CASE
	call lu(model, model, 'L' )			!L: lower case, U: upper case

	Ising= .FALSE.;  XYZ = .FALSE.
	if(model.eq.'ising') Ising = .TRUE.
	if(model.eq.'xyz') XYZ = .TRUE.

	read(0, *) tmcs					! total Monte Carlo Steps per spins
	read(0, *) tmcs_eq				! total Monte Carlo Steps per spins for equilibration
	read(0, *) ht, lt, tint				! higher temp, lower temp., temp. interval
	read(0, *) (s(i), i = 1, nspecies)		! spin-moment
	read(0, *) n_speci_incl				! no. of species to include

	allocate(species_to_include(n_speci_incl), stgg(n_speci_incl))

	read(0, *) (species_to_include(i), i = 1, n_speci_incl) ! species lable to include in calculation
	read(0, *) sc(1), sc(2), sc(3)			! supercell dimension

	if(mod(sc(1), 2).ne.0.or.mod(sc(2), 2).ne.0.or.mod(sc(3), 2).ne.0) then
		write(6, *) ''
		write(6, *) '	ERROR'
		write(6, *) '	~~~~~'
		write(6, *) ''
		write(6, *) '	Dimension of supercell along x, y, and z must be in the multiples of 2'
		write(6, *) '	STOPPING now'
		write(6, *) ''
		stop
	end if

	read(0, *) staggered				! logic for staggered magnetization
	read(0, *) bc(1), bc(2), bc(3)			! boundary condition along x, y, z

	! CONVERTING STRING VARIBLES INTO LOWER CASE
	call lu(bc(1), bc(1), 'L' )
	call lu(bc(2), bc(2), 'L' )
	call lu(bc(3), bc(3), 'L' )

	read(0, *) sample 				! sample per unit MC calculations)

	read(0, *) angle, dphi	 			! logic, least angle window (helps in convergence)

	dphi = convert_to_rad*dphi

	read(0, *) to_cal 				! at this step interval all observables will be calculated
	read(0, *) Zeeman, h(1:3)			! Magnetic field logic, Mx, My, Mz (in T)
	read(0, *) g_factor				! g-factor
	read(0, *) single_ion_anisotropy		! Logic(.T./.F.)
	read(0, *) para, para_value			! Logic, parameter value --> calculation in terms of parameter
        read(0, *) ovrr, overr_para, ovrr_start_interval ! overrelaxed (T/F), overrelaxed parameter (0-1), interval of overrelaxation starts
        overr_para = abs(overr_para)
        read(0, *) seed_value				! Seed for random number generation
	seed = seed_value
        read(0, *) dope					! Logic for dopping

	!read(0, *) EXalgo, exchange_interval, temp_ex	!Exchange algo (logic), interval to exchange replicas, for using given equi-spaced temprature (T/F)
	!read(0, *) beta_critria, optbeta, beta_file	!Beta(M) convergence criteria, optimization steps for beta(M), if Beta.dat file exists(logic)

	if((para_value.eq.0).and.para.and.rank == 0) then
		write(6, *) ""
		write(6, *) "ERROR:: parameter value is set ON and value &
       			cannot be zero. Kindly have a look into the input file"
       		write(6, *) ""
		STOP
	end if

	close(0)					! closing of input file
	
        end subroutine read_input
