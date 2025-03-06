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

	subroutine get_nbd_extents
	
	use init

	implicit none
	
	integer :: i, j, nbd_dis_max

	logical :: file_present
		
	inquire(file='j_exchange', exist=file_present)
	if(.not.file_present) then
	if (rank == 0) then
		write(6, *) "==> Required file 'j_exchange' is not present"
		write(6, *) "	 STOPPING now"
		write(6, *) ""
	end if
		stop
	end if

	open(2, file='j_exchange', status='old', action='read')
	read(2, *) no_of_nbd					! no. of diff. bond length w.r.t. any central ion
	if(allocated(nbd_dis)) deallocate(nbd_dis)
	allocate(nbd_dis(no_of_nbd))				! nbd distances

	do i = 1, no_of_nbd
		read(2, *) nbd_dis(i), similar_bonds
		do j = 1, similar_bonds
			read(2, *)
		end do
	end do

	!closing of j_exchange file
	close(2)

	nbd_dis_max = maxval(nbd_dis, dim=1)

	! Extent of nbd cells around a lattice point
	nbd_cell_x = ceiling(nbd_dis_max/lp(1))
	nbd_cell_y = ceiling(nbd_dis_max/lp(2))
	nbd_cell_z = ceiling(nbd_dis_max/lp(3))

	! Simulation box start from_ till to_
	fromx = nbd_cell_x + 1; tox = sc(1) + nbd_cell_x
	fromy = nbd_cell_y + 1; toy = sc(2) + nbd_cell_y
	fromz = nbd_cell_z + 1; toz = sc(3) + nbd_cell_z

	end subroutine get_nbd_extents
