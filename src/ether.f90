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

	! DATE & TIME
	character(8) :: date
	character(10) :: time
	character(5) :: zone

	integer, dimension(8) :: start, finish
	integer :: days, hrs, mins, secs, n_seed, &
		stepi, ei, ej, el
	integer, dimension(:), allocatable :: seed

	! Initialize MPI
	call MPI_INIT(ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

	call random_seed(size=n_seed)
	allocate(seed(n_seed))

	! initialize the random number seed = 1992
	seed = 1992
    
	call date_and_time(VALUES=start)

	if (rank == 0) call system('rm -f *.dat *.xsf')

	call get_tot_species(nspecies, total_ions_per_cell)

	allocate(s(nspecies), sia_vec(1:3, nspecies), &
		x(total_ions_per_cell), y(total_ions_per_cell), z(total_ions_per_cell), &
		species(0:nspecies), ionn(0:nspecies))

	call read_input
	call random_seed(put=seed)
	if (rank == 0) call startup(start)
	call read_structure(nspecies, total_ions_per_cell, lp, abc, x, y, z, ionn, species, rank)
	call j_values
	call get_sia_values
	call parameters
	call generate_supercell
	if (rank == 0) call write_initial_conf
	call getting_nbd
	if (rank == 0) call write_nbd
	call generate_bubble_indices

	! Total temperature count
	nscan = nint((ht - lt) / (tint)) + 1

	call allocate_observables

	! For initiating Ground Spin State (GSS) file
	if (rank == 0) call write_gss(0)
	! For creating output files
	if (rank == 0) call write_output_files(0)

	allocate(addmoreitr(size), num_iterations(size), istart(0:size), iend(0:size), &
		tlobs(size), tlspn(size))

	addmoreitr = 0; istart = 0; iend = 0
     
	! Calculate the reminder (left) values is size divides the nscan,
	! it can be used for assigning this leftovers to other processors.
	left = mod(nscan, size)
	do el = 1, left

		addmoreitr(el) = 1  ! add more iterations

	end do

	! Calculate the work distribution: number of iterations per process
	interval = nscan / size
	do ei = 1, size

		num_iterations(ei) = addmoreitr(ei) + interval

	end do

	! Calculate the starting and end points of distributed temperature index
	do ei = 1, size

		istart(ei) = iend(ei-1) + 1
		iend(ei) = iend(ei-1) + num_iterations(ei)

	end do
    
	! Number of observable (like: temp., energy, Cv, Mag., Chi,..., etc.)
	total_observables = 14

	! Allocate array for each process
	! local observable length
    	local_olen = total_observables + 3 * nspecies  		! total observable + 
                                                     		! 3 magnetic component * no. of species
    	local_slen = product(sc) * lattice_per_unit_cell * 4  	! total no. of cell * lattice per unit cel * 
                                                          	! (3 spin vectors + ID)

	do ei = 1, size

		! Total length of observables for local array
		tlobs(ei) = local_olen*(iend(ei) - istart(ei) + 1)

		! Total length of spin details for local array
		tlspn(ei) = local_slen*(iend(ei) - istart(ei) + 1)

	end do
	
	allocate(local_obs(tlobs(rank+1)), local_spn(tlspn(rank+1)), si_obs(size), &
		si_spn(size))

	local_obs = 0d0; local_spn = 0d0
    	
	! Calculate the starting index for recieving buffer from non-root processors
	si_obs(1) = 0
	si_spn(1) = 0
	do ei = 2, size

		! For observables
		si_obs(ei) = sum(tlobs(1:ei-1))
		! For spin details
		si_spn(ei) = sum(tlspn(1:ei-1))
		
        end do

    	if (rank == 0) then

        !$OMP PARALLEL
        num_of_threads = omp_get_num_threads()
        !$OMP END PARALLEL

        write(6, *) ''
        write(6,'(" Calculations are running on total: ")')
        write(6, *) ''
        write(6,'("              OpenMP threads =>", i3)') num_of_threads
        write(6,'("              MPI processes  =>", i3)') size
        write(6, *) ''

	write(6, "(' List of temperatures')")
	write(6, "(' ~~~~~~~~~~~~~~~~~~~~')")
	write(6, *) ''

	do ei = 0, size-1

		allocate(temp_assigned(num_iterations(ei + 1)))
		el = 0
		do ej = istart(ei + 1), iend(ei + 1)
			el = el + 1
			temp_assigned(el) = ht - tint * (ej - 1)
		end do
		write(6, "(' At process ID', i3, ' are:')") ei
		write(6, "('                      ',5g11.4)") temp_assigned
		write(6, *) ''
		deallocate(temp_assigned)

	end do
		
	! Allocate global array only on the root process (rank 0)
    	allocate(global_obs(local_olen*nscan), global_spn(local_slen*nscan))

	else
		! This is necessay because for command line flag -check=all it gives error
		! so we are making these buffer as dummy things for non-root processes
		allocate(global_obs(0), global_spn(0))
	end if

	el = 0
	li_obs = 0	! local index for observables
	li_spn = 0	! local index for spins
	temperature: do itemp = istart(rank + 1), iend(rank + 1)

		! local index for local_* array
		el = el + 1
		li_obs = local_olen*(el - 1)
		li_spn = local_slen*(el - 1)

		! Temperature (T)
		temp = ht - tint * (itemp - 1)
		! kBT^-1
		beta = 1.0d0 / (kb*temp)

		! Creating *spK* file
		initiate_spin_files = .TRUE.
		call write_spins_at_K(itemp)

		! Setting zero to sample variables
		call zeroes('sample')

		acceptance_counting = 0
		sampling: do samplei = 1, sample

			call fresh_spins
			call zeroes('eng_mag')

			do stepi = 1, tmcs

			if (ovrr.and..not.Ising.and. &
				(mod(real(stepi), real(ovrr_start_interval)) == 0.0)) then

				call overrelaxation

			end if

			call mc(stepi, acceptance_count)
			acceptance_counting = acceptance_counting + acceptance_count

			end do

			call evaluate_observables('avg_eng')
			call evaluate_observables('avg_mag')

		end do sampling

        	call evaluate_observables('observables_per_sample')

        	! At temperature K

        	! Getting magnetic moment vectors
        	call get_moment_vectors
        	! Evaluate all observables
        	call evaluate_observables('store')

        	! Writing spin states into *spK* file
        	call write_spins_at_K(itemp)

		! Writing Job status
		write(6, "(' Task for temperature ', g11.4, 'is completed.')") temp

	end do temperature

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

		call evaluate_observables('write')
		call graph

		call system('rm -rf data spins')
		call system('mkdir spins')
		call system('mkdir data')
		call system('mv *spK* spins')
		call system('mv *conf.xsf *.dat data')
		call system('mv graph* data')

		write(6, *) ''
		write(6, *) "May the force be with you"
		write(6, *) "~ Mukesh Kumar Sharma"
		write(6, *) "e-mail@ msharma1@ph.iitr.ac.in"
		write(6, *) ''
		write(6, *) ''

		close(10004)
		close(10005)

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

