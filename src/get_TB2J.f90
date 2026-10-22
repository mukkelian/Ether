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

	subroutine get_TB2J()
	
	use init, only: dp, nspecies, J_iso, J_prime, J_ani, DMI, &
				no_of_nbd, nbd_dis, species, Jiso, Jprime, DMvec, Jani, &
				nbd_finding_criteria
	use mpi
	implicit none

	real(dp) :: vec(3)
	integer :: i, j, k, l, m, n, J_ID(2), ios
	logical :: exists, file_exists(4)
	character (len=20) :: files(4)
	character (len=3) :: atom_i, atom_j, atom(2), out
	character (len=80) :: line
	files(1) = 'J_iso'
	files(2) = 'J_ani'
	files(3) = 'J_prime'
	files(4) = 'DMI'

	! For TB2J
	nbd_finding_criteria = 0.001
	file_exists = .False.
	do i = 1, 4
		exists = .False.
		inquire(file=trim(files(i)), exist=exists)
		file_exists(i) = exists
	end do
	J_iso = file_exists(1)
	J_ani = file_exists(2)
	J_prime = file_exists(3)
	DMI = file_exists(4)

	if (nspecies==0) call terminate('Number of species should not be zero')
		! reading J_iso file
	if(J_iso) then
		open(20, file= 'J_iso', status='old', action='read')
		do i = 1, 3;read(20, *);end do
		read(20, '(a)') line
		read(line, *) no_of_nbd
	end if

		! reading J_ani file
	if(J_ani) then
		open(21, file= 'J_ani', status='old', action='read')
		do i = 1, 3;read(21, *);end do
		read(21, '(a)') line
		read(line, *) no_of_nbd
	end if

	! reading J_prime file
	if(J_prime) then
		open(22, file= 'J_prime', status='old', action='read')
		do i = 1, 3;read(22, *);end do
		read(22, '(a)') line
		read(line, *) no_of_nbd
	end if

	! reading DMI file
	if(DMI) then
		open(23, file= 'DMI', status='old', action='read')
		do i = 1, 3;read(23, *);end do
		read(23, '(a)') line
		read(line, *) no_of_nbd
	end if

	if (no_of_nbd==0) call terminate('Bonds are absent')
	allocate(nbd_dis(no_of_nbd))
	if (J_iso) allocate (Jiso(no_of_nbd, nspecies, nspecies, -1:1, -1:1, -1:1))
	if (J_ani) allocate (Jani(1:3, 1:3, no_of_nbd, nspecies, nspecies, -1:1, -1:1, -1:1))
	if (J_prime) allocate (Jprime(no_of_nbd, nspecies, nspecies, -1:1, -1:1, -1:1))
	if (DMI) allocate (DMvec(1:3, no_of_nbd, nspecies, nspecies, -1:1, -1:1, -1:1))

	n = no_of_nbd
	! reading J_iso file
	if(J_iso) then
		Jiso = 0.0_dp
		do i = 1, n;read(20, *) nbd_dis(i); end do

		do
		read(20, '(a)', iostat=ios, end=10) line
		read(line, *) atom_i, atom_j, i, vec(1:3), l, m, n
		call process_label(atom_i); call process_label(atom_j)

		atom = (/ atom_i, atom_j /)
		do j = 1, 2
			j_ID(j) = 0
			do k = 1, nspecies
				if (atom(j).eq.species(k)) then
					j_ID(j) = k
				end if
			end do
			if (j_ID(j) == 0) call terminate('species not found')
		end do
		read(20, *) Jiso(i, j_ID(1), J_ID(2), l, m, n)
		end do
10	close(20)
	end if

	! reading J_ani file
	if(J_ani) then
		Jani = 0.0_dp
		do i = 1, n;read(21, *) nbd_dis(i); end do

		do
			read(21, '(a)') line
			read(line, *) atom_i, atom_j, i, vec(1:3), l, m, n
			call process_label(atom_i); call process_label(atom_j)

			atom = (/ atom_i, atom_j /)
			do j = 1, 2
				do k = 1, nspecies
					if (atom(j).eq.species(k)) then
						j_ID(j) = k
					end if
				end do
			end do
			read(21, *) Jani(1, 1:3, i, j_ID(1), J_ID(2), l, m, n)
			read(21, *) Jani(2, 1:3, i, j_ID(1), J_ID(2), l, m, n)
			read(21, *, end=11) Jani(3, 1:3, i, j_ID(1), J_ID(2), l, m, n)
		end do
11	continue
	close(21)
	end if

	! reading J_prime file
	if(J_prime) then
		Jprime = 0.0_dp
		do i = 1, n;read(22, *) nbd_dis(i); end do

		do
			read(22, '(a)') line
			read(line, *) atom_i, atom_j, i, vec(1:3), l, m, n
			call process_label(atom_i); call process_label(atom_j)

			atom = (/ atom_i, atom_j /)
			do j = 1, 2
				do k = 1, nspecies
					if (atom(j).eq.species(k)) then
						j_ID(j) = k
					end if
				end do
			end do
			read(22, *, end=12) Jprime(i, j_ID(1), J_ID(2), l, m, n)
		end do
12	continue
	close(22)
	end if

	! reading DMI file
	if(DMI) then
		DMvec = 0.0_dp
		do i = 1, n;read(2, *) nbd_dis(i); end do

		do
			read(23, '(a)') line
			read(line, *) atom_i, atom_j, i, vec(1:3), l, m, n
			call process_label(atom_i); call process_label(atom_j)

			atom = (/ atom_i, atom_j /)
			do j = 1, 2
				do k = 1, nspecies
					if (atom(j).eq.species(k)) then
						j_ID(j) = k
					end if
				end do
			end do
			read(23, *, end=13) DMvec(1:3, i, j_ID(1), J_ID(2), l, m, n)
		end do
13	continue
	close(23)
	end if

	contains

	subroutine process_label(lbl)

	implicit none
	character(len=*), intent(inout) :: lbl

	call lu(lbl(1:1), out, 'U' )
	lbl(1:1) = out
	call lu(lbl(2:2), out, 'L' )
	lbl(2:2) = out

	end subroutine process_label

	end subroutine get_TB2J
