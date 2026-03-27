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

	subroutine write_temperature_infos
	
	use init
	use omp_lib

	implicit none

	integer :: ei, ej, el
 
        !$OMP PARALLEL
        num_of_threads = omp_get_num_threads()
        !$OMP END PARALLEL

        write(6, *) ''
        write(6,'(" Calculations are running on total: ")')
        write(6, *) ''
        write(6,'("              OpenMP threads =>", i3)') num_of_threads
        write(6,'("              MPI processes  =>", i3)') nprocs
        write(6, *) ''

	write(6, "(' List of temperatures')")
	write(6, "(' ~~~~~~~~~~~~~~~~~~~~')")
	write(6, *) ''

	do ei = 0, nprocs-1

		allocate(temp_assigned(int(total_temperatures(ei + 1))))
		el = 0
		do ej = istart(ei + 1), nscan, nprocs
			el = el + 1
			temp_assigned(el) = temperature(ej)
		end do
		write(6, "(' At process ID', i3, ' are:')") ei
		write(6, "(24X,5g11.4)") temp_assigned
		write(6, *) ''
		deallocate(temp_assigned)

	end do
	
	end subroutine write_temperature_infos
