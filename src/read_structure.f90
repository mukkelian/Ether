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

	subroutine read_structure(nsp, nions, lp, abc, coordx, coordy, coordz, ions, species_name, rank_)

		implicit none
		
		integer :: tions, i, j
		integer, intent(in) :: nsp, nions, rank_
		integer, intent(out) :: ions(0:nsp)

		real(8), intent(out) :: lp(3), abc(3, 3), coordx(nions), coordy(nions), coordz(nions)
		real(8) :: wt
		
		character(len=200) :: title, label
		character(len=2), intent(out) :: species_name(0:nsp)
		
		logical :: file_present

		inquire(file='structure.vasp', exist=file_present)
		if(.not.file_present) then
			write(6, *) "==> 'structure.vasp' is not present"
			write(6, *) "	  STOPPING now"
			write(6, *) ""
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
		ions = 0
		read(10001,*) (ions(i), i =1, nsp)
		tions = sum(ions)

		coordx = 0; coordy = 0; coordz = 0
		read(10001, *) label
		if(trim(adjustl(label)).eq.'Direct')then
			write(6, *) 'Choose POSCAR in cartesian co-ordinate only'
			write(6, *) 'STOPPING'
			stop
		endif

		do i = 1, tions
			read(10001,*) coordx(i), coordy(i), coordz(i) !cartesian co-ordinates
		end do
		close(10001)
	    	if (rank_ == 0) write(6, *) '==> Reading structure is completed!'
		
	end subroutine read_structure
