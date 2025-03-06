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

	subroutine overrelaxation

	use init

	implicit none

	integer :: i, j, k, l, m

	real(dp) :: A_ovr(3), Si(3), rn_ovr

        overrelaxation_method : do m = 1, overr_steps
        
        	call get_random_indices(i, j, k, l)

		A_ovr = real(0.0, dp)

		call get_ovrr_vec(i, j, k, l, A_ovr)
	

		Si(1:3) = ion(1:3, i, j, k, l)

	        Si(1:3) = (2* &
	        (dot_product(A_ovr(1:3), Si(1:3))/dot_product(A_ovr(1:3), A_ovr(1:3)) )* &
	        A_ovr(1:3)) - Si(1:3)

		ion(1:3, i, j, k, l) = Si(1:3)/sqrt(dot_product(Si(1:3), Si(1:3)))

	end do overrelaxation_method

	end subroutine overrelaxation
