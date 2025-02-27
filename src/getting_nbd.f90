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

	integer :: search_nbd_x(2), search_nbd_y(2), search_nbd_z(2), i, j, k, l, &
		ii, jj, kk, ll, nbd_count, nbd_dis_max, central_ion_ID, nbd_ID, ith_bond
	
	real(8) :: distance, origin(3), nbd_pos(3)

	allocate(nn(no_of_nbd, total_ions, 0:20, 0:4))

	if (rank == 0) write(6, *) '==> Sensing neighbourhoods'
	nn = 0

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(sc, lattice_per_unit_cell, no_of_nbd, ion, nn, nbd_dis, &
!$OMP& nbd_cell_x, nbd_cell_y, nbd_cell_z, nbd_finding_criteria)
	!$OMP DO SCHEDULE(dynamic)
	central_atom :	do k = 1, sc(3) + 2*nbd_cell_z
			call nbd_finding_limit(sc(3), nbd_cell_z, k, search_nbd_z)

			do j = 1, sc(2) + 2*nbd_cell_y
			call nbd_finding_limit(sc(2), nbd_cell_y, j, search_nbd_y)

			do i = 1, sc(1) + 2*nbd_cell_x
			call nbd_finding_limit(sc(1), nbd_cell_x, i, search_nbd_x)

			do l = 1, lattice_per_unit_cell

			central_ion_ID = int(ion(0, i, j, k, l))
			origin = ion(6:8, i, j, k, l)

			! no. of distinct bond
			nbd : do ith_bond = 1, no_of_nbd
			nbd_count = 0.0
							
			do kk = search_nbd_z(1), search_nbd_z(2)
				do jj = search_nbd_y(1), search_nbd_y(2)
					do ii = search_nbd_x(1), search_nbd_x(2)
						do ll = 1, lattice_per_unit_cell

			nbd_pos = ion(6:8, ii, jj, kk, ll)
			nbd_ID = int(ion(0, ii, jj, kk, ll))

			if(central_ion_ID.ne.nbd_ID)then                 

			distance = sqrt(dot_product(origin - nbd_pos, origin - nbd_pos ))
										
				if(abs(nbd_dis(ith_bond) - (distance)).le.nbd_finding_criteria)then

					nbd_count = nbd_count +1
					! Storing nbd's ID				
					nn(ith_bond, central_ion_ID, nbd_count, 0) = nbd_ID

					! Storing nbd's indices
					nn(ith_bond, central_ion_ID, nbd_count, 1) = ii
					nn(ith_bond, central_ion_ID, nbd_count, 2) = jj
					nn(ith_bond, central_ion_ID, nbd_count, 3) = kk

					! Storing nbd's atom count
					nn(ith_bond, central_ion_ID, nbd_count, 4) = ll
			
				end if        
               		end if

						end do
					end do
				end do
			end do

			! Storing total no. of similar
			! nbds on particular ion for mth distinct bond
			nn(ith_bond, central_ion_ID, 0, 0) = nbd_count

			end do nbd								
			end do
	end do
	end do
	end do central_atom
        !$OMP END DO
!$OMP END PARALLEL
	if (rank == 0) write(6, *) '==> Neighbourhood informations are created now'

	contains
	
	subroutine nbd_finding_limit(total_extension, nbd_upto, cell_pos, to_nbd_cell)

		implicit none

		integer, intent (in) :: total_extension, nbd_upto, cell_pos
		integer, intent (out) :: to_nbd_cell(2)

			to_nbd_cell(1) = cell_pos - nbd_upto
			to_nbd_cell(2) = cell_pos + nbd_upto
			if(cell_pos .le. nbd_upto) to_nbd_cell(1) = 1
			if(cell_pos .gt. (total_extension + 2*nbd_upto - nbd_upto)) &
			to_nbd_cell(2) = total_extension + 2*nbd_upto
			
	end subroutine nbd_finding_limit
	
	end subroutine getting_nbd
