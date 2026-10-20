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

!=======================================================================
!  vector_direction.f90
!
!  Gives the direction (sign) of every component of a vector:
!
!      vec = ( 2.0304,  3.44, -2.334 )   ->   dir = ( 1, 1, -1 )
!
!  Arguments
!  ---------
!   n    (in)  : integer, number of components (3 for a 3D vector)
!   vec  (in)  : real(8) array of size n, the vector
!   dir  (out) : integer array of size n, the direction (-1, 0 or 1)
!=======================================================================
	subroutine vector_direction(n, vec, dir)
	implicit none

	integer, parameter :: dp = kind(1.0d0)

	integer,  intent(in)  :: n
	real(dp), intent(in)  :: vec(n)
	integer,  intent(out) :: dir(n)

	real(dp), parameter :: tol = 0.0_dp
	integer :: i

	do i = 1, n
		if (vec(i) > tol) then
			dir(i) = 1
		else if (vec(i) < -tol) then
			dir(i) = -1
		else
			dir(i) = 0
		end if
	end do

	end subroutine vector_direction
