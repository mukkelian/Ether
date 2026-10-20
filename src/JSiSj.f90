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

	!Hamiltonian: JSiSj
	subroutine JSiSj(on, i, Si, central_ion_ID, total_energy)
		use init, only: dp, ion, nn, j_exc, no_of_nbd
		implicit none

		integer, intent(in) :: i, central_ion_ID
		integer :: total_nbr, ID_num, ion_ID, cell, inbr, nbdi

		real(dp), intent(in) :: Si(3)
		real(dp) :: Sj(3), SiSj(3), Jij(3), eout
		real(dp), intent(inout) :: total_energy

		logical, intent(in) :: on	! to calculate the energy on i-th lattice site

		ion_ID = int(ion(0, i))
		! Over distinct bonds
		do nbdi = 1, no_of_nbd

			! For nbdi bond total connecting neighbours to the central ion [ID = ion(0, i)] is
			total_nbr = nn(nbdi, ion_ID, 0, 0)

			! no. of similar nbd for ith distinct bond (nbdi)
			do inbr = 1, total_nbr

				! nbr's ID in the cell
				cell = nn(nbdi, ion_ID, inbr, 1)
				ID_num = int(ion(4, cell))

				if(ID_num.eq.0) then
					go to 4
				end if

				Sj(1:3) = ion(1:3, cell)

				!Si.Sj
				SiSj = Si*Sj	!point-wise multiplication

				!Jij term
				Jij = j_exc(nbdi, central_ion_ID, ID_num, 1:3)

				if(on) then
					!Jij.Si.Sj term, for on site energy, no double counting
					eout = dot_product(SiSj, Jij)
				else
					!(1/2)*Jij.Si.Sj term, double counting is considered in calculating total energy
					eout = 0.5_dp*dot_product(SiSj, Jij)
				end if

				total_energy = total_energy + eout

4				continue
			end do
		end do

	end subroutine JSiSj
