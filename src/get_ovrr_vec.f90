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

	subroutine get_ovrr_vec(io, ovrr_vec)

	use init

	implicit none

	integer, intent(in) :: io
	integer :: ionID

	real(dp), intent(inout) :: ovrr_vec(3)

	ovrr_vec = 0.0_dp

	! Species ID of central ion
	ionID = int(ion(4, io))
	
	call JSiSj_ovrr(io, ionID, ovrr_vec)

	end subroutine get_ovrr_vec

	! JSiSj	
	subroutine JSiSj_ovrr(i, central_ion_ID, ovrr_vec)

	use init, only: dp, nn, ion, j_exc, no_of_nbd

	implicit none

	integer, intent(in) :: i, central_ion_ID
	integer :: total_nbr, ion_ID, cell, ith_bond, ith_nbr

	real(dp) :: Sj(3), Jij(3), JSj(3)
	real(dp), intent(inout) :: ovrr_vec(3)

	ion_ID = int(ion(0, i))

	! Over distinct bonds
	do ith_bond = 1, no_of_nbd

		! For ith_bond bond total connecting neighbours to the central ion [ID = ion(0, i)] is
		total_nbr = nn(ith_bond, ion_ID, 0, 0)

		! no. of similar nbd for ith distinct bond (ith_bond)
		do ith_nbr = 1, total_nbr

			! nbr's ID in the cell
			cell = nn(ith_bond, ion_ID, ith_nbr, 1)

			if(ion(4, cell).eq.0) then
				go to 4
			end if

			Sj(1:3) = ion(1:3, cell)

			!Jij term
			Jij = j_exc(ith_bond, central_ion_ID, int(ion(4, cell)), 1:3)

			!Jij.Sj term
			JSj = Jij*Sj

			ovrr_vec = ovrr_vec + JSj

4			continue
		end do
	end do

	end subroutine JSiSj_ovrr
