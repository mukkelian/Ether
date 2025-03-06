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

	subroutine George_Marsaglia(spinx, spiny, spinz)

		use init, only: dp

		implicit none

		real(dp), intent(out) :: spinx, spiny, spinz
		real(dp) :: v1, v2, S_val, rn

7		call get_random_num(-real(1, dp), real(1, dp), rn)
		v1 = rn
		call get_random_num(-real(1, dp), real(1, dp), rn)
		v2 = rn
		S_val = v1**2 + v2**2
		
		if (S_val.ge.1) goto 7
		
		spinx = 2*v1*sqrt(1-S_val)
		spiny = 2*v2*sqrt(1-S_val)
		spinz = 1-2*S_val

	end subroutine George_Marsaglia
