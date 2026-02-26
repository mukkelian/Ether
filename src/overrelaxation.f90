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

	integer :: i, m

	real(dp) :: A_ovr(3), Si(3), rn_ovr, on_site_eng

        overrelaxation_method : do m = 1, ovrr_steps
        
        	call get_random_indices(total_ions, i)

		A_ovr = 0.0_dp

		call get_ovrr_vec(i, A_ovr)

		Si(1:3) = ion(1:3, i)

		! check energy before applying overrelaxation algorithm
		!call Hamiltonian(.FALSE., i, Si, Si, on_site_eng)
		!print*, 'BEFORE OVRR:'
		!print*, 'Eng:',on_site_eng
		!print*, 'Central ION vec:', Si(1:3)

	        Si(1:3) = (2* &
	        (dot_product(A_ovr(1:3), Si(1:3))/dot_product(A_ovr(1:3), A_ovr(1:3)) )* &
	        A_ovr(1:3)) - Si(1:3)

		ion(1:3, i) = Si(1:3)/sqrt(dot_product(Si(1:3), Si(1:3)))

		! check energy before applying overrelaxation algorithm
		!call Hamiltonian(.FALSE., i, ion(1:3, i), ion(1:3, i), on_site_eng)
		!print*, 'AFTER OVRR:'
		!print*, ' Eng:', on_site_eng
		!print*, 'Central ION vec:', ion(1:3, i)
	
	end do overrelaxation_method

	end subroutine overrelaxation
