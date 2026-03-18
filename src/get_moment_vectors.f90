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

	subroutine get_moment_vectors

	use init, only: dp, ion, s, total_ions, mm_vector

	implicit none
		
	integer :: i, j
		
	mm_vector = 0.0_dp

	do i = 1, total_ions
                j = int(ion(4, i))
	        ! moment vectors for each magentic ion
	        mm_vector(j, 1:3) = &
		        mm_vector(j, 1:3) &
		        + ion(1:3, i)*s(j)
	end do

	end subroutine get_moment_vectors
