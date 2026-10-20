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

	!Hamiltonian: gmbSH	
	subroutine gmbSH(Si, central_Ion, total_energy)
		use init, only: dp, mb, g_factor, H, s, scaled, &
			ScaledSpin, xyz
		implicit none

		integer, intent(in) :: central_Ion
		real(dp) :: eout
		real(dp), intent(in) :: Si(3)
		real(dp), intent(inout) :: total_energy

		! energy due to magnetic field
		if(scaled.and.xyz) then
			eout = -g_factor*mb*ScaledSpin(central_Ion)*dot_product(Si, H)
		else
			eout = -g_factor*mb*s(central_Ion)*dot_product(Si, H)
		end if
		total_energy = total_energy + eout
	end subroutine gmbSH
