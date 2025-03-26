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

	subroutine get_tot_species(n, total_ions_in_cell)
	
		implicit none
		
		integer :: i, j, k, jj
		integer, intent(out) :: n, total_ions_in_cell
		integer, allocatable :: totion(:)

		character(len=200) :: titel

		logical :: file_found

		inquire(file='structure.vasp', exist=file_found)
		if(.not.file_found) then
			write(6, *) "==> 'structure.vasp' is not present"
			write(6, *) "	  STOPPING now"
			write(6, *) ""
			stop
		end if

		open(unit=0, file='structure.vasp', status='old', action='read')

        n = 0

		do jj = 1, 6
			read(0, *) 	! skipping 6 lines 
		end do

        read(0, '(a)') titel

        call count_species(titel, n)

		allocate(totion(n))
		read(titel, *) totion	! total ions
		total_ions_in_cell = sum(totion)

		close(0)

	end subroutine get_tot_species
