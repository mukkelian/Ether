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

	subroutine get_vector_direction

	use init, only: dp, ion, no_of_nbd, nn, nbd_capacity, total_ions, vd
	implicit none

	integer :: i
	integer :: total_nbr, ID_num, ion_ID, cell, inbr, nbdi, dir(3)

	real(dp) :: rj(3), ri(3), rij(3)


	if(allocated(vd)) deallocate(vd)
	if(total_ions == 0) call terminate ('Zero magnetic ions')
	allocate(vd(1:3, nbd_capacity, total_ions))
	vd = 0
	do i = 1, total_ions
	ion_ID = int(ion(0, i))
	ri = ion(6:8, i)
	! Over distinct bonds
	do nbdi = 1, no_of_nbd
		! For nbdi bond total connecting neighbours to the central ion [ID = ion(0, i)] is
		total_nbr = nn(nbdi, ion_ID, 0, 0)

		! no. of similar nbd for ith distinct bond (nbdi)
		do inbr = 1, total_nbr

			! nbr's ID in the cell
			cell = nn(nbdi, ion_ID, inbr, 1)
			ID_num = int(ion(4, cell))

			if(ID_num.eq.0) then
				go to 4
			end if
			! Rj
			rj(1:3) = ion(6:8, cell)

			!Rij vector
			rij = rj - ri
			call vector_direction(3, rij, dir)
			vd(1:3, nbdi, i) = dir(1:3)

4			continue
		end do
	end do
	end do

	end subroutine get_vector_direction
