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

	subroutine search_nbd(pos0, pos1, nbd_range, captured)

	use init

	implicit none

	real(dp), intent (in), dimension(3) :: pos0, pos1
	real(dp), intent(in) :: nbd_range
	real(dp) :: R_vec(3), distance, ion_pos(3)
	integer :: shift(3), i, j, k, ii, jj, kk, &
		to_x, to_y, to_z
	logical, intent (out) :: captured

	captured = .FALSE.
	to_x = 3; to_y = 3; to_z = 3 
	
	if(bc(1).eq.'o') to_x = 1
	if(bc(2).eq.'o') to_y = 1
	if(bc(3).eq.'o') to_z = 1

	shift = (/0, 1, -1/)
	do i = 1, to_x
	ii = shift(i)
		do j = 1, to_y
		jj = shift(j)
			do k = 1, to_z
			kk = shift(k)

	ion_pos = SCabc(1, 1:3)*ii &
		+ SCabc(2, 1:3)*jj &
		+ SCabc(3, 1:3)*kk &
		+ pos1(1:3)

	R_vec = pos0 - ion_pos
	distance = sqrt(dot_product(R_vec, R_vec))
	if(abs(nbd_range - distance).le.nbd_finding_criteria) then
		captured = .TRUE.
		return
	end if

			end do
		end do
	end do

	end subroutine search_nbd
