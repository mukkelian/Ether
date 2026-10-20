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

	!Hamiltonian: siaS2
	subroutine siaS2(Si2, central_Ion, total_energy)
		use init, only: dp, sia_vec
		implicit none

		integer, intent(in) :: central_Ion
		real(dp) :: sia(3)
		real(dp), intent(in) :: Si2(3)
		real(dp), intent(inout) :: total_energy

		sia(1:3) = sia_vec(1:3, central_Ion)

		! energy due to single ion anisiatropy
		total_energy = total_energy + dot_product(sia, Si2)

	end subroutine siaS2
