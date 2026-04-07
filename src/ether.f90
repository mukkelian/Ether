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
	integer :: comm_ei, color, dim1(2), dim2(2)
	
	logical :: file_found
	
	character(len=10) :: Date, Month

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
	deallocate (x, y, z)
	if(ssp) call get_spiral_state(0.0_dp, 'ions_for_spiral')

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

	! Calculate the starting and end points of distributed temperature index
	allocate(total_temperatures(nprocs), istart(0:nprocs))
	do ei = 1, nprocs
		istart(ei) = ei
		el = 0
		do ej = istart(ei), nscan, nprocs
			el = el + 1
		end do
		total_temperatures(ei) = el 
	end do
    
	! Number of observable (like: temp., energy, Cv, Mag., Chi,..., etc.)
	total_observables = 20

    	if (rank == 0) call write_temperature_infos

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

		! Temperature (T)
		temp = temperature(itemp)
		! kBT^-1
		beta = 1.0_dp/(kb*temp)

		! Creating *spK* file
		initiate_spin_files = .TRUE.
		call write_spins_at_K(temp, itemp)

		! Setting zero to repeat variables
		call zeroes('repeat')

		acceptance_counting = 0

		! SAMPLING
		sampling: do repeati = 1, repeat

                        call random_seed(put=seed+itemp+175*rank+666*repeati)
			call fresh_spins
			call zeroes('avg')

		swap_count = 0
		perform_MCS: do stepi = 1, tmcs

		total_energy = 0.0_dp

		if (ovrr.and.XYZ.and.&
			.not.Zeeman.and.&
			.not.SIA.and.&
			(mod(real(stepi), real(ovrr_MCS)) .eq. 0).and.&
			(mod(real(stepi), real(to_cal)) .ne. 0)) then
			call overrelaxation
		end if

		! Monte Carlo with Metropolis aLgo
		call Metropolis(beta, acceptance_count)
		acceptance_counting = acceptance_count + &
				acceptance_counting

		! Calculate Total Energy
		call get_tot_energy(total_energy)

		! For getting the status of MCS simulations.
    		inquire(file='STATUS', exist=file_found)
    		if (file_found) call get_status(rank, stepi, repeati, temp)
    		
		! Wait for other MPIs to complete
		call MPI_Barrier(comm_ei, ierr)

		if (file_found .and. rank.eq.0) call system('rm -f STATUS')

		! Monte Carlo with Parallel Tempering (PT) algo
		PT: if(PTalgo) then
			if (mod(stepi, exchange_interval) == 0) then
				swap_count = swap_count + 1

    			call parallel_tempering(ion, total_energy, beta, itemp, &
        			swap_count, comm_ei)

			end if 
		end if PT

		! Calculate observable (Total Energy, Magnetisation, etc.) after Equilibration
		call observable(beta, stepi)

		end do perform_MCS

		call get_eng(beta, 'avg_eng')
		call get_mag(beta, 'avg_mag')
		if(ssp) call get_spiral_state(beta, 'avg_spiral_state')

		end do sampling

		! Acceptance_counting per repeat*tmcs
		acceptance_counting = acceptance_counting/(repeat*tmcs)

        	call observables_per_repeat

        	! At temperature K

        	! Getting magnetic moment vectors
        	call get_moment_vectors

                ! Store all calculated observables
        	call process_observables('store', itemp)

        	! Writing spin states into *spK* file
        	call write_spins_at_K(temp, itemp)
        	
        	! Writing spin states into ETHER.spn
        	dim1 = (/0, total_info/)
        	dim2 = (/1, total_ions/)
        	call rw_file('w', 'ETHER.spn', 502, dim1, &
        		dim2, itemp, ion)

		! Writing Job status
		write(6, "(' Task for temperature '&
			, g11.4, 'is completed.')") temp
    		call MPI_Comm_free(comm_ei, ierr)
    		
1		continue

	end do temperature_MPI

	deallocate(ion)

	write(6, "(' Process ID', i3, ' is free now, assigned work is completed.')") rank

	! Root process will finalize observables and write results after all data has been received
	if (rank == 0) then

		! Write all collected observables
		call process_observables('write', 0)
		call graph

		call system('mv *spK* spins')
		call system('mv *network.xsf *.dat data')
		call system('mv graph* data')

		write(6, *) ''
		write(6, *) ''
		write(6, *) 'SUMMARY :::::::::::::::::::::::::::::::::::::'
		write(6, *) ''
		write(Date, '(i0)') start(3)
		write(Month, '(i0)') start(2)
		write(6, '("        PROGRAM STARTED on date ",A,"-",A,"-",i4)') &
			trim(Date), trim(Month), start(1)
		write(6, "('        at time ',i2,' hrs. ',i2,' min. ',i2,' sec. ')") &
			start(5:7)
			
		call date_and_time(VALUES=finish)

		write(6, *) ''
		write(Date, '(i0)') finish(3)
		write(Month, '(i0)') finish(2)
		write(6, '("        PROGRAM ENDED on date ",A,"-",A,"-",i4)') &
			trim(Date), trim(Month), finish(1)
		write(6, "('        at time ',i2,' hrs. ',i2,' min. ',i2,' sec. ')") &
			finish(5:7)
		write(6, *) ''
		write(6, *) ':::::::::::::::::::::::::::::::::::::::::::::'
		write(6, *) ''

	end if

	! Finalize MPI
	call MPI_FINALIZE(ierr)

	end program ether

