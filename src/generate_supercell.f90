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

	subroutine generate_supercell
	
	use init
	
	implicit none

	integer :: i, j, k, l, s_count

	real(dp), allocatable :: atom(:, :), ar(:, :)
	real(dp) :: sx, sy, sz, rnum
	logical :: file_found

	! SAVING DETAILS OF ATOMS PRESENT IN STRUCTURE FILE
	allocate(atom(total_ions_per_cell, 0:8))
	atom = 0
	atomic_details: do k = 1, nspecies
		do i = int(sum(ionn(0:k-1)))+1, int(sum(ionn(0:k)))

			atom(i,0) = i			! atom no.
			if(Ising) then
				atom(i, 3) = (-1)**(i)
			else
				if (angle.eqv..FALSE.) then
					call George_Marsaglia(sx, sy, sz)
					atom(i, 3) = sz	! S(z)
					atom(i, 1) = sx	! S(x)
					atom(i, 2) = sy	! S(y)
					atom(i, 5) = 0	! phi value
				else

					call get_random_num(-real(1, dp), real(1, dp), rnum)
					atom(i, 3) = rnum						! S(z)
					call get_random_num(real(0, dp), 2*pi, rnum)
          				atom(i, 1) = sqrt(real(1, dp)-(atom(i, 3)**2))*cos(rnum)	! S(x)
					atom(i, 2) = sqrt(real(1, dp)-(atom(i, 3)**2))*sin(rnum)	! S(y)
					atom(i, 5) = rnum					! phi value
				end if
			endif
			atom(i, 4) = k		! tag for element
			atom(i, 6) = x(i)	! x
			atom(i, 7) = y(i)	! y
			atom(i, 8) = z(i)	! z

 		enddo
	enddo atomic_details

	! FILTRATION OF SPECIES FOR CALC.
	s_count = 0
	do i = 1, n_speci_incl
		do j =1, nspecies
			if (species_to_include(i).eq.species(j)) then
				s_count = s_count + ionn(j)
			end if
		end do
	end do

	allocate(ar(s_count, 0:8), stgg_ion(s_count))
	ar = 0

	! STAGGERED INFO
	if (staggered) call get_staggered_info(s_count)

	! Getting total lattice per unit cell, which will be used for MC calculations
	lattice_per_unit_cell = 0
	inclusion : do i = 1, n_speci_incl
		do j =1, nspecies
			
		if (species_to_include(i).eq.species(j)) then
		do k = sum(ionn(0:j-1))+1, sum(ionn(0:j))
			lattice_per_unit_cell = lattice_per_unit_cell +1	! Counting included ions
			do l = 0, 8
			ar(lattice_per_unit_cell, l) = atom(k, l)		! 'k' will be used for identifying
			end do 							! these element in further last moment
		enddo
		endif

		enddo
	enddo inclusion

	deallocate(atom)

	! EXPANDING INTO SUPERCELL
	allocate(ion(0:8, sc(1) + 2*nbd_cell_x, sc(2) + 2*nbd_cell_y, &
		sc(3) + 2*nbd_cell_z, lattice_per_unit_cell))
	s_count = 0; ion = 0

	do l = 1, lattice_per_unit_cell
		do k = 1, sc(3) + 2*nbd_cell_z
			do j = 1, sc(2) + 2*nbd_cell_y
				do i = 1, sc(1) + 2*nbd_cell_x
			


		s_count = s_count +1				! Counting total no. of ions in supercell
		ion(0, i, j, k, l) = s_count			! ID
		ion(1:5, i, j, k, l) = ar(l, 1:5)		! Sx, Sy, Sz, species no. & phi
		ion(6:8, i, j, k, l) = abc(1, 1:3)*(i - 1) &	! co-ordinates
					+ abc(2, 1:3)*(j - 1) + abc(3, 1:3)*(k - 1) &
					+ ar(l, 6:8)

				end do
			end do
		end do
	end do

	deallocate(ar)

	! total no. of lattice points
	total_ions = s_count	! in super cell + in boundary layer
	total_lattice_sites = product(sc)*lattice_per_unit_cell	! in super cell == simulation box

	if (ovrr) overr_steps = nint(total_lattice_sites*overr_para)
	if (rank == 0) then
		write(6, *) '==> Supercells are generated'
		if (ovrr) write(6, "(' ==> Total number of OVRR steps: ',I6)") overr_steps
		if (ovrr) write(6, "('     MCS steps after which OVRR will be performed: ',I4)") ovrr_start_interval
	end if
	end subroutine generate_supercell

