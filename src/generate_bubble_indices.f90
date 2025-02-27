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

	subroutine generate_bubble_indices

		use init
	    
		implicit none
	    
		allocate(bblx(sc(1)), bbly(sc(2)), bblz(sc(3)))

		bblx = 0; bbly = 0; bblz = 0

		call get_index(bblx, nbd_cell_x + 1, fromx, tox)
		call get_index(bbly, nbd_cell_y + 1, fromy, toy)
		call get_index(bblz, nbd_cell_z + 1, fromz, toz)

	contains
	
	subroutine get_index(indices, stepi, from, to)
	
		implicit none

		integer :: i, j, k, extension
		integer, intent(in) :: stepi, from, to
		integer, allocatable, intent(out) :: indices(:)

		extension = to - from + 1
		allocate(indices(extension))
		j = from; k = 0
		do i = 1, extension

			indices(i) = j + stepi*k
			 k = k + 1
			if(indices(i).gt.to) then
				j = j + 1; k = 0
				indices(i) = j
				k = k + 1
			end if

		end do

	end subroutine get_index
	
	end subroutine generate_bubble_indices
