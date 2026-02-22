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

	subroutine getting_nbd

	use init
	use omp_lib

	implicit none

	integer :: search_nbd_x(2), search_nbd_y(2), search_nbd_z(2), i, j, &
		nbd_count, nbd_dis_max, central_ion_ID, nbd_ID, ith_bond
	
	real(dp) :: distance, origin(3), nbd_pos(3), R_vec(3)

	integer :: nbd_capacity
	
	logical :: captured = .FALSE.

	nbd_capacity = 20
	allocate(nn(no_of_nbd, total_ions, 0:nbd_capacity, 0:1))

	if (rank == 0) write(6, *) '==> Sensing neighbourhoods'
	nn = 0

	SCabc(1, :) =  sc(1)*abc(1, :)
	SCabc(2, :) =  sc(2)*abc(2, :)
	SCabc(3, :) =  sc(3)*abc(3, :)

	!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(sc, lattice_per_unit_cell, no_of_nbd, ion, nn, nbd_dis, &
	!$OMP& total_ions, nbd_finding_criteria, SCabc)
	!$OMP DO SCHEDULE(dynamic)
	central_atom :	do i = 1, total_ions

		central_ion_ID = int(ion(0, i))
		origin = ion(6:8, i)

		! no. of distinct bond
		nbd : do ith_bond = 1, no_of_nbd
		nbd_count = 0

		! nbd ion
		do j = 1, total_ions
		nbd_pos = ion(6:8, j)
		nbd_ID = int(ion(0, j))

		if(central_ion_ID.ne.nbd_ID)then
			call search_nbd(origin, nbd_pos, captured)						
			if(captured)then

				nbd_count = nbd_count + 1
				! Storing nbd's ID				
				nn(ith_bond, central_ion_ID, nbd_count, 0) = nbd_ID

				! Storing nbd's indices
				nn(ith_bond, central_ion_ID, nbd_count, 1) = j
			
			end if        
       		end if

		end do

		! Storing total no. of similar
		! nbds on particular ion for mth distinct bond
		nn(ith_bond, central_ion_ID, 0, 0) = nbd_count

		end do nbd

	end do central_atom
        !$OMP END DO
	!$OMP END PARALLEL

	if (rank == 0) write(6, *) '==> Neighbourhood informations are created now'
	
	end subroutine getting_nbd
