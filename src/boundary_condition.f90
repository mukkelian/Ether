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

	subroutine boundary_condition(ith, jth, kth, lth)

		use init

		implicit none

		integer, intent(in) :: ith, jth, kth, lth

		! BOUNDARY X
		along_x : if(bc(1).eq.'c')then
			if (ith.eq.sc(1)+1) then
				ion(1:5, 1, jth, kth, lth) = &
				ion(1:5, ith, jth, kth, lth)
			elseif(ith.eq.2)then
				ion(1:5, sc(1)+2, jth, kth, lth) = &
				ion(1:5, ith, jth, kth, lth)
			end if
		end if along_x
                                                                
		! BOUNDARY Y
		along_y : if(bc(2).eq.'c')then
			if (jth.eq.sc(2)+1) then
				ion(1:5, ith, 1, kth, lth) = &
				ion(1:5, ith, jth, kth, lth)
			elseif(jth.eq.2)then
				ion(1:5, ith, sc(2)+2, kth, lth) = &
				ion(1:5, ith, jth, kth, lth)
			end if
		end if along_y

                 ! BOUNDARY Z
                  along_z : if(bc(3).eq.'c')then
			if (kth.eq.sc(3)+1) then
				ion(1:5, ith, jth, 1, lth) = &
				ion(1:5, ith, jth, kth, lth)
			elseif(kth.eq.2)then
				ion(1:5, ith, jth, sc(3)+2, lth) = &
				ion(1:5, ith, jth, kth, lth)
			end if
		end if along_z

	end subroutine boundary_condition
