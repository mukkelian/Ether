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

	subroutine get_ovrr_vec(io, jo, ko, lo, ovrr_vec)

		use init

		implicit none

		integer, intent(in) :: io, jo, ko, lo
		integer :: ionID, nbdi, inbr, cell

		real(8) :: Si_ref(3), sia_value(3)
		real(8), intent(inout) :: ovrr_vec(3)

		ovrr_vec = 0

		! Central ion's species ID
		ionID = int(ion(4, io, jo, ko, lo))
		
		! Central ion's spin vectors
		Si_ref = ion(1:3, io, jo, ko, lo)

		call JSiSj_ovrr(io, jo, ko, lo, ionID, ovrr_vec)

		if(Zeeman) call gmbSH_ovrr(g_factor, mb, H, ovrr_vec)

		if(single_ion_anisotropy) then
			sia_value(1:3) = sia_vec(1:3, ionID)
			call siaS2_ovrr(Si_ref, sia_value, ovrr_vec)
		end if

	contains

	! JSiSj	
	subroutine JSiSj_ovrr(i, j, k, l, central_ion_ID, ovrr_vec)
	
		implicit none

		integer, intent(in) :: i, j, k, l, central_ion_ID
		integer :: total_nbr, posx, posy, posz

		real(8) :: Sj(3), Jij(3), JSj(3)
		real(8), intent(inout) :: ovrr_vec(3)

		! Over distinct bonds
		do nbdi = 1, no_of_nbd

			! For nbdi bond total connecting neighbours to the central ion [ID = ion(0, i, j, k, l)] is
			total_nbr = nn(nbdi, int(ion(0, i, j, k, l)), 0, 0)

			! no. of similar nbd for ith distinct bond (nbdi)
			do inbr = 1, total_nbr

				! nbr's cell position
				posx = nn(nbdi, int(ion(0, i, j, k, l)), inbr, 1)
				posy = nn(nbdi, int(ion(0, i, j, k, l)), inbr, 2)
				posz = nn(nbdi, int(ion(0, i, j, k, l)), inbr, 3)

				! nbr's ID in the cell
				cell = nn(nbdi, int(ion(0, i, j, k, l)), inbr, 4)

				if(ion(4, posx, posy, posz, cell).eq.0) then
					go to 4
				end if

				Sj(1:3) = ion(1:3, posx, posy, posz, cell)

				!Jij term
				Jij = j_exc(nbdi, central_ion_ID, int(ion(4, posx, posy, posy, cell)), 1:3)

				!Jij.Sj term
				JSj = Jij*Sj

				ovrr_vec = ovrr_vec + JSj

4				continue
			end do
		end do

	end subroutine JSiSj_ovrr
	
	!gmbSH	
	subroutine gmbSH_ovrr(g, mu, H, ovrr_vec)

		implicit none

		real(8) :: ref_vec(3)
		real(8), intent(in) :: g, mu, H(3)
		real(8), intent(inout) :: ovrr_vec(3)

		! OVRR vec. due to magnetic field
		ref_vec = -(g*mu*H)
		
		ovrr_vec = ovrr_vec + ref_vec
	
	end subroutine gmbSH_ovrr

	!SINGLE ION ANISOTRPY (SIA)
	subroutine siaS2_ovrr(Si, sia_val, ovrr_vec)

		implicit none

		real(8) :: ref_vec(3)
		real(8), intent(in) :: Si(3), sia_val(3)
		real(8), intent(inout) :: ovrr_vec(3)

		! OVVR vec. due to single ion anisiatropy
		ref_vec = Si*sia_val			

		ovrr_vec = ovrr_vec + ref_vec

	end subroutine siaS2_ovrr

	end subroutine get_ovrr_vec
