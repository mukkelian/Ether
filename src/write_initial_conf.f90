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

		integer :: i, j, k, l, sce(3)

		open(unit=2, file='starting_spin_conf.xsf', status='unknown')

		write(2, *) 'CRYSTAL'
		write(2, *) 'PRIMVEC'
		sce = sc	! Supercell extent (sce)
		if(bc(1).eq.'o') sce(1) = sc(1)+1
		if(bc(2).eq.'o') sce(2) = sc(2)+1
		if(bc(3).eq.'o') sce(3) = sc(3)+1
		write(2, *) (abc(1,j)*(sce(1)), j= 1,3)
		write(2, *) (abc(2,j)*(sce(2)), j= 1,3)
		write(2, *) (abc(3,j)*(sce(3)), j= 1,3)
		write(2, *) 'CONVEC'
		write(2, *) (abc(1,j)*(sce(1)), j= 1,3)
		write(2, *) (abc(2,j)*(sce(2)), j= 1,3)
		write(2, *) (abc(3,j)*(sce(3)), j= 1,3)
		write(2, *) 'PRIMCOORD'
		write(2, *) product(sc)*lattice_per_unit_cell, " 1"

		! SPINS
		do l = 1, lattice_per_unit_cell
			do k = fromz, toz
				do j = fromy, toy
					do i = fromx, tox


		if (ion(4, i, j, k, l).ne.0) then

			write(2, 101) species(int(ion(4, i, j, k, l))), &
			ion(6:8, i, j, k, l), ion(1:3, i, j, k, l)

		end if

					end do
				end do
			end do
		end do
101	   	format(A5,7f13.7)

		close(2)

		write(6, *) "==> Initial spin states have be written into 'starting_spin_conf.xsf'"

	end subroutine write_initial_conf
