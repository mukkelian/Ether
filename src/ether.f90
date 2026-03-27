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

	program ether

	use init
	use omp_lib
	use mpi

	implicit none

	integer, dimension(8) :: start, finish
	integer :: n_seed, stepi, ei, ej, el, max_temp_count, swap_count
	integer :: comm_ei, color
	real(dp) :: beta_before
	! Initialize MPI
	call MPI_INIT(ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

	call random_seed(size=n_seed)
	allocate(seed(n_seed))

	! initialize the random number with seed = 1992
	seed = 1992
    
	call date_and_time(VALUES=start)

	if (rank == 0) call system('rm -f *.dat *.xsf')

	call get_tot_species(nspecies, total_ions_per_cell)

	allocate(s(nspecies), &
		x(total_ions_per_cell), &
		y(total_ions_per_cell), &
		z(total_ions_per_cell), &
		species(0:nspecies), ions(0:nspecies))

	call read_input
	call random_seed(put=seed+rank)
	if (rank == 0) call startup(start)
	call read_structure(nspecies, total_ions_per_cell, lp, abc, &
		x, y, z, ions, species, rank)
	call generate_supercell
	call j_values
	call get_sia_values
	call parameters
	if (rank == 0) call write_fetched_lattice_network
	call getting_nbd
	if (rank == 0) call write_nbd

	! Total temperature count
	nscan = nint((ht - lt) / (tint)) + 1

        allocate(temperature(nscan))
        do ei = 1, nscan
        	temperature(ei) =  ht - tint * (ei - 1)
        end do

	call allocate_observables

	! For initiating Ground Spin State (GSS) file
	if (rank == 0) call write_gss(0)
	! For creating output files
	if (rank == 0) call write_output_files(0)

	allocate(addmoreitr(nprocs), num_iterations(nprocs), istart(0:nprocs), &
		tlobs(nprocs), tlspn(nprocs))

	addmoreitr = 0; istart = 0

	! Calculate the starting and end points of distributed temperature index
	allocate(total_temperatures(nprocs))
	do ei = 1, nprocs
		istart(ei) = ei
		el = 0
		do ej = istart(ei), nscan, nprocs
			el = el + 1
		end do
		total_temperatures(ei) = el 
	end do
    
	! Number of observable (like: temp., energy, Cv, Mag., Chi,..., etc.)
	total_observables = 14

	! Allocate array for each process
	! local observable length
    	local_olen = total_observables + 3*nspecies  		! total observable + 
                                                     		! 3 magnetic component * no. of species
    	local_slen = product(sc)*lattice_per_unit_cell*4  	! total no. of cell * lattice per unit cel * 
                                                          	! (3 spin vectors + ID)

	do ei = 1, nprocs

		! Total length of observables for local array
		tlobs(ei) = local_olen*total_temperatures(ei)

		! Total length of spin details for local array
		tlspn(ei) = local_slen*total_temperatures(ei)

	end do
	
	allocate(local_obs(tlobs(rank+1)), local_spn(tlspn(rank+1)), si_obs(nprocs), &
		si_spn(nprocs))

	local_obs = 0.0_dp; local_spn = 0.0_dp
    	
	! Calculate the starting index for recieving buffer from non-root processors
	si_obs(1) = 0
	si_spn(1) = 0
	do ei = 2, nprocs
		! For observables
		si_obs(ei) = sum(tlobs(1:ei-1))
		! For spin details
		si_spn(ei) = sum(tlspn(1:ei-1))
        end do

    	if (rank == 0) then
		call write_temperature_infos
		! Allocate global array only on the root process (rank 0)
    		allocate(global_obs(local_olen*nscan), global_spn(local_slen*nscan))
	else
		! This is necessay because for command line flag -check=all it gives error
		! so we are making these buffer as dummy things for non-root processes
		allocate(global_obs(0), global_spn(0))
	end if

	li_obs = 0	! local index for observables
	li_spn = 0	! local index for spins

	max_temp_count = maxval(total_temperatures)
	allocate(temp_range(nprocs, max_temp_count))
	temp_range = 0
	do ei = 1, nprocs
		el = 0
		do ej = istart(ei), nscan, nprocs
			el = el + 1
			temp_range(ei, el) = ej
		end do
	end do

	el = 0

	temperature_MPI: do ei = 1, max_temp_count
	
		itemp = temp_range(rank + 1, ei)

		! Create dynamic communicator for this iteration
		if (itemp /= 0) then
			color = ei
		else
			color = MPI_UNDEFINED
		end if

		call MPI_Comm_split(MPI_COMM_WORLD, color, rank, comm_ei, ierr)

		if(itemp.eq.0) goto 1

		! local index for local_olen/local_slen array
		el = el + 1
		li_obs = local_olen*(el - 1)
		li_spn = local_slen*(el - 1)

		! Temperature (T)
		temp = temperature(itemp)
		! kBT^-1
		beta = 1.0_dp/(kb*temp)

		! Creating *spK* file
		initiate_spin_files = .TRUE.
		call write_spins_at_K(itemp)

		! Setting zero to repeat variables
		call zeroes('repeat')

		acceptance_counting = 0

		! SAMPLING
		sampling: do repeati = 1, repeat

                        call random_seed(put=seed+itemp+175*rank+666*repeati)
			call fresh_spins
			call zeroes('eng_mag')

		swap_count = 0
		perform_MCS: do stepi = 1, tmcs

		total_energy = 0.0_dp

		if (ovrr.and.XYZ.and.&
			.not.Zeeman.and.&
			.not.SIA.and.&
			(mod(stepi, ovrr_MCS) == 0)) then
			call overrelaxation
		end if

		! Monte Carlo with Metropolis aLgo
		call Metropolis(beta, acceptance_count)
		acceptance_counting = acceptance_count + &
				acceptance_counting

		! Calculate Total Energy
		call get_tot_energy(total_energy)

		! Wait for other MPIs to complete
		call MPI_Barrier(comm_ei, ierr)

		! Monte Carlo with Parallel Tempering (PT) algo
		PT: if(PTalgo) then
			if (mod(stepi, exchange_interval) == 0) then
				swap_count = swap_count + 1

    			call parallel_tempering(ei, ion, total_energy, beta, itemp, &
        			swap_count, comm_ei)

			end if 
		end if PT

		! Calculate observable (Total Energy, Magnetisation) after Equilibration
		call evaluate_observable(beta, stepi)

		end do perform_MCS

		call evaluate_eng_observables(beta, 'avg_eng')
		call evaluate_mag_observables(beta, 'avg_mag')

		end do sampling

		! Acceptance_counting per repeat*tmcs
		acceptance_counting = acceptance_counting/(repeat*tmcs)

        	call evaluate_observables_per_repeat

        	! At temperature K

        	! Getting magnetic moment vectors
        	call get_moment_vectors

                ! Evaluate all observables
        	call process_observables('store')

        	! Writing spin states into *spK* file
        	call write_spins_at_K(itemp)

		! Writing Job status
		write(6, "(' Task for temperature '&
			, g11.4, 'is completed.')") temp

    		call MPI_Comm_free(comm_ei, ierr)

1		continue
	end do temperature_MPI

	! Gather data from all processes to root
	call MPI_GATHERV(local_obs, tlobs(rank+1), MPI_DOUBLE_PRECISION, &
			global_obs, tlobs, si_obs, &
			MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

	call MPI_GATHERV(local_spn, tlspn(rank+1), MPI_DOUBLE_PRECISION, &
			global_spn, tlspn, si_spn, &
			MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

	write(6, "(' Process ID', i3, ' is free now, assigned work is completed.')") rank

	! Root process will finalize observables and write results after all data has been received
	if (rank == 0) then

		call process_observables('write')
		call graph

		call system('rm -rf data spins')
		call system('mkdir spins')
		call system('mkdir data')
		call system('mv *spK* spins')
		call system('mv *network.xsf *.dat data')
		call system('mv graph* data')

		write(6, *) ''
		write(6, *) ''
		write(6, *) 'SUMMARY :::::::::::::::::::::::::::::::::::::'
		write(6, *) ''
		write(6, '("        PROGRAM STARTED on date ",i2,"-",i2,"-",i4)') start(3), start(2), start(1)
		write(6, "('        at time ',i2,' hrs. ',i2,' min. ',i2,' sec. ')") start(5:7)
		call date_and_time(VALUES=finish)
		write(6, *) ''
		write(6, '("        PROGRAM ENDED on date ",i2,"-",i2,"-",i4)') finish(3), finish(2), finish(1)
		write(6, "('        at time ',i2,' hrs. ',i2,' min. ',i2,' sec. ')") finish(5:7)
		write(6, *) ''
		write(6, *) ':::::::::::::::::::::::::::::::::::::::::::::'
		write(6, *) ''

	end if

	! Finalize MPI
	call MPI_FINALIZE(ierr)

	end program ether

