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

		use init

		implicit none
		
		integer :: i, j, k, l
		
		mm_vector = 0

		do i = fromx, tox
			do j = fromy, toy
				do k = fromz, toz
					do l = 1, lattice_per_unit_cell
						! moment vectors for each magentic ion
						mm_vector(int(ion(4, i, j, k, l)), 1:3) = &
						mm_vector(int(ion(4, i, j, k, l)), 1:3) &
						+ ion(1:3, i, j, k, l)*s(int(ion(4, i, j, k, l)))
					end do
				end do
			end do
		end do

	end subroutine get_moment_vectors
