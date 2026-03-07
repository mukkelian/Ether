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

	subroutine write_nbd

		use init

		implicit none
	
		integer :: i, m

		open(unit=2, file='nbd.dat', status='unknown')
		write(6, *) "==> Neighbourhood informations are writing into the 'nbd.dat' file"
                write(2, *) "# This file can be used to confirm the neighbourhood details of any selected"
                write(2, *) "# central ion with ion number (<ION no.>)."
                write(2, *) "# [NOTE: The content of nbd.dat is independent of the species labels present in 'structure.vasp'."
                write(2, *) "# The contents of nbd.dat depend only on the bond lengths provided in the 'j_exchange' file.]"
                write(2, *) "# To visualize the structure, use the 'VESTA' software (https://jp-minerals.org/vesta/en/download.html)"
                write(2, *) "# and open the generated 'fetched_lattice_network.xsf' file."
                write(2, *) "# Check the neighbourhood details (given in this file) using the <Atom Label> for the selected central <ION no.> IDs."
		write(2, *) ''

		nbd_inf : do i = 1, total_ions

			write(2,10) int(ion(0, i))
			write(2,*)"~~~~~~~~~~~~~~"
 
			nbd_write : do m = 1, no_of_nbd
				write(2, 11) nbd_dis(m), &
				nn(m, int(ion(0, i)), &
				1:nn(m, int(ion(0, i)), 0, 0), &
				 0)

10	      	format(" ION no. ",i5)
11      	format(" For bond length: ",f12.8,", nbds are ==> "20i7)
			end do nbd_write

			write(2,*)''

		end do nbd_inf

		close(2)

	end subroutine write_nbd
