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

        subroutine cross_product(vec1, vec2, cp)

        use init, only: dp
        
        implicit none

        real(dp), dimension(3), intent(in) :: vec1, vec2
	real(dp), dimension(3), intent(out) :: cp
	
		cp(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
		cp(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
		cp(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

	end subroutine cross_product
