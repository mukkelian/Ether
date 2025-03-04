! Code Ether, based on the Monte Carlo technique, can be used to
! study the static and dynamics of spin models applied to any ion geometry.
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
		if (ith.gt.sc(1)) then
			ion(1:5, ith - sc(1), jth, kth, lth) = &
			ion(1:5, ith, jth, kth, lth)
		elseif(ith.lt.(fromx + nbd_cell_x))then
			ion(1:5, ith + sc(1), jth, kth, lth) = &
			ion(1:5, ith, jth, kth, lth)
		end if
	end if along_x

	! BOUNDARY Y
	along_y : if(bc(2).eq.'c')then
		if (jth.gt.sc(2)) then
			ion(1:5, ith, jth - sc(2), kth, lth) = &
			ion(1:5, ith, jth, kth, lth)
		elseif(jth.lt.(fromy + nbd_cell_y))then
			ion(1:5, ith, jth + sc(2), kth, lth) = &
			ion(1:5, ith, jth, kth, lth)
		end if
	end if along_y

	! BOUNDARY Z
	along_z : if(bc(3).eq.'c')then
		if (kth.gt.sc(3)) then
			ion(1:5, ith, jth, kth - sc(3), lth) = &
			ion(1:5, ith, jth, kth, lth)
		elseif(kth.lt.(fromz + nbd_cell_z))then
			ion(1:5, ith, jth, kth + sc(3), lth) = &
			ion(1:5, ith, jth, kth, lth)
		end if
	end if along_z

	end subroutine boundary_condition
