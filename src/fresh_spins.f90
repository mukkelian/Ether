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

	subroutine fresh_spins

		use init
	
		implicit none
	
		integer :: i
		real(dp) :: S_vec_previous(5), S_vec_updated(5)

		S_vec_previous = 0; S_vec_updated = 0

		do i = 1, total_ions

			S_vec_previous = ion(1:5, i)

			call update_spin_details(S_vec_previous, S_vec_updated)
			ion(1:5, i) = S_vec_updated

		end do

	end subroutine fresh_spins
