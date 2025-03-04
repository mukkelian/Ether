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

	subroutine update_spin_details(ref1, ref2)

		use init

		implicit none

		real(8) :: rnum, sx, sy, sz, phi
		real(8), intent(in) :: ref1(5)
		real(8), intent(out) :: ref2(5)
		
		ref2 = ref1							! storing species no.
		if(Ising) then
			ref2(3) = -ref1(3)
		else
			if(angle.eqv..FALSE.) then
				call George_Marsaglia(sx, sy, sz)
				ref2(1) = sx; ref2(2) = sy; ref2(3) = sz
				ref2(5) = 0
			else
				ref2(5) = ref1(5)				! old phi
				call get_random_num(-1d0, +1d0, rnum)
				ref2(3) = rnum					! cos(theta) : S(z)
				call get_random_num(-1d0, +1d0, rnum)
				phi = ref2(5) + (dphi*rnum)			! older_phi + delta*rnum
				ref2(1) = sqrt(1d0-(ref2(3)**2))*cos(phi)	! S(x)
				ref2(2) = sqrt(1d0-(ref2(3)**2))*sin(phi)	! S(y)
				ref2(5) = phi					! storing phi value
			end if
		end if

		! Normalising to mod 1
		ref2(1:3) = ref2(1:3)/sqrt(dot_product(ref2(1:3), ref2(1:3)))

	end subroutine update_spin_details
