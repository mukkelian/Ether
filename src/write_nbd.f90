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
	
		integer :: i, j, k, l, m

		open(unit=2, file='nbd.dat', status='unknown')
		write(6, *) "==> Neighbourhood informations are writing into the 'nbd.dat' file"

		write(2, *) "This file can be used to confirm the neighbourhood details of any central ion."
		write(2, *) "To do so, use 'XCrySDen' (http://www.xcrysden.org/) and open the generated"
		write(2, *) "'starting_spin_conf.xsf' as file > Open Structure > Open XSF > starting_spin_conf.dat."
		write(2, *) "Check nbd's details by clicking the <Atom Info> (at bottom) and see the 'Selected Atom No.' informations. "
		write(2, *) ''

		nbd_inf : do i = 1, sc(1) + 2*nbd_cell_x
				do j = 1, sc(2) + 2*nbd_cell_y
					do k = 1, sc(3) + 2*nbd_cell_z
	                        		do l = 1, lattice_per_unit_cell

						write(2,10) int(ion(0, i, j, k, l))
						write(2,*)"~~~~~~~~~~~~~~"
 
						nbd_write : do m = 1, no_of_nbd
							write(2, 11) m, &
							nn(m, int(ion(0, i, j, k, l)), &
							1:nn(m, int(ion(0, i, j, k, l)), 0, 0), &
							 0)

10	      	format(" ION no. ",i5)
11      	format(" For bond length no. ",i2," nbds are ==> "20i7)
						end do nbd_write

						write(2,*)''

                                                end do
					end do
				end do
		end do nbd_inf

		close(2)

	end subroutine write_nbd
