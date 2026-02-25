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

	subroutine read_structure(nsp, nions, lp, abc, &
		coordx, coordy, coordz, ions, species_name, rank)

	use init, only: dp

	implicit none

	integer :: tions, i, j
	integer, intent(in) :: nsp, nions, rank
	integer, intent(out) :: ions(0:nsp)

	real(dp), intent(out) :: lp(3), abc(3, 3), coordx(nions), &
		coordy(nions), coordz(nions)
	real(dp) :: wt, frac_pos(3), cart_pos(3)
		
	character(len=200) :: title, label
	character(len=2), intent(out) :: species_name(0:nsp)
	character(len=2) :: ab
	
	logical :: file_present

	inquire(file='structure.vasp', exist=file_present)
	if(.not.file_present) then
	if (rank == 0) then
		write(6, *) "==> Required file 'structure.vasp' is not present"
		write(6, *) "	 STOPPING now"
		write(6, *) ""
	end if
		stop
	end if

	open(unit=10001, file='structure.vasp', status='old', action='read')
               read(10001,*) title
        read(10001,*) wt

	lattice_constant : do i = 1,3
		read(10001,*) (abc(i,j), j = 1,3)
	enddo lattice_constant

	lp(1) = sqrt(dot_product(abc(1, 1:3), abc(1, 1:3)))	! lattice parameter a
	lp(2) = sqrt(dot_product(abc(2, 1:3), abc(2, 1:3)))	! lattice parameter b
	lp(3) = sqrt(dot_product(abc(3, 1:3), abc(3, 1:3)))	! lattice parameter c

	read(10001,*) (species_name(i), i = 1,nsp)
	species_name(0) = 'X'

	do i = 1, nsp
		ab = species_name(i)
		call lu(ab(1:1), ab(1:1), 'U' )
		call lu(ab(2:2), ab(2:2), 'L' )
		species_name(i) = ab
	end do

	ions = 0
	read(10001,*) (ions(i), i =1, nsp)
	tions = sum(ions)

	coordx = 0; coordy = 0; coordz = 0
	read(10001, *) label
	
	call lu(label, label, "L")

	if(trim(adjustl(label)).eq.'direct')then

		do i = 1, tions
			! fractional co-ordinates
			read(10001,*) coordx(i), coordy(i), coordz(i)
			frac_pos = (/coordx(i), coordy(i), coordz(i)/)
			call frac_to_cart(abc, frac_pos, cart_pos)
			coordx(i) = cart_pos(1)
			coordy(i) = cart_pos(2)
			coordz(i) = cart_pos(3)
		end do
 		close(10001)

	elseif(trim(adjustl(label)).eq.'cartesian')then
	
		! cartesian co-ordinates
		do i = 1, tions
			read(10001,*) coordx(i), coordy(i), coordz(i)
		end do
		close(10001)

	else
		if (rank == 0) print*, ''
		if (rank == 0) print*, '    Incorrect structure file'
		if (rank == 0) print*, ''
		if (rank == 0) print*, '    STOPPING now'
		stop
	end if

    	if (rank == 0) write(6, *) '==> Reading structure is completed!'
		
	end subroutine read_structure
