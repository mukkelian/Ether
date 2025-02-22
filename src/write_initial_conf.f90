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

	! WRITING INITIAL STRUCTURE FILE
	subroutine write_initial_conf

		use init

		implicit none

		integer :: i, j, k, l	
	
		open(unit=10003, file='starting_spin_conf.xsf', status='unknown')
		write(10003, *) 'ATOM'
		do i = 1, sc(2) + 2*nbd_cell_x
			do j = 1, sc(2) + 2*nbd_cell_y
				do k = 1, sc(3) + 2*nbd_cell_z
					do l = 1, lattice_per_unit_cell

						write(10003,10021) species(int(ion(4, i, j, k, l))), &
	                                        ion(6:8, i, j, k, l), ion(1:3, i, j, k, l)

					end do
				end do
			end do
		end do
		close(10003)

10031		format(12A7)
10021   	format(A5,7f13.7)

		write(6, *) "==> Initial spin states have be written into 'starting_spin_conf.dat'"

	end subroutine write_initial_conf
