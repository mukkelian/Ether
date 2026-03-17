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
	subroutine write_fetched_lattice_network

		use init

		implicit none

		integer :: i, j

		open(unit=2, file='fetched_lattice_network.xsf', status='unknown')

		write(2, *) 'CRYSTAL'
		write(2, *) 'PRIMVEC'
		write(2, *) (abc(1,j)*(sc(1)), j= 1,3)
		write(2, *) (abc(2,j)*(sc(2)), j= 1,3)
		write(2, *) (abc(3,j)*(sc(3)), j= 1,3)
		write(2, *) 'CONVEC'
		write(2, *) (abc(1,j)*(sc(1)), j= 1,3)
		write(2, *) (abc(2,j)*(sc(2)), j= 1,3)
		write(2, *) (abc(3,j)*(sc(3)), j= 1,3)
		write(2, *) 'PRIMCOORD'
		write(2, *) product(sc)*lattice_per_unit_cell, " 1"

		! SPINS
		do i = 1, total_ions

		if (ion(4, i).ne.0) then
			write(2, 101) species(int(ion(4, i))), &
			ion(6:8, i)!, ion(1:3, i)
		end if

		end do
101	   	format(A5,7f13.7)

		close(2)

		write(6, *) "==> Fetched lattice networks have be written into 'fetched_lattice_network.xsf'"

	end subroutine write_fetched_lattice_network
