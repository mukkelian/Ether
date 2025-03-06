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

	subroutine Hamiltonian(ih, jh, kh, lh, Si_center, total_ene)
	
		use init

		implicit none

		integer, intent(in) :: ih, jh, kh, lh
		integer :: ionID, nbdi, inbr, cell
		
		real(dp) :: sia_value(3)
		real(dp), intent(in) :: Si_center(3)	! central ion
		real(dp), intent(out) :: total_ene

		total_ene = 0

		! Central ion's species ID
		ionID = int(ion(4, ih, jh, kh, lh))

		call JSiSj(ih, jh, kh, lh, Si_center, ionID, total_ene)

		if(Zeeman) call gmbSH(Si_center, g_factor, mb, H, total_ene)

		if(single_ion_anisotropy.and..not.Ising) then
			sia_value(1:3) = sia_vec(1:3, ionID)
			call siaS2(Si_center, sia_value, total_ene)
		end if
	
	contains

	! JSiSj	
	subroutine JSiSj(i, j, k, l, Si, central_ion_ID, total_energy)
	
		implicit none

		integer, intent(in) :: i, j, k, l, central_ion_ID
		integer :: total_nbr, posx, posy, posz, ID_num, ion_ID

		real(dp), intent(in) :: Si(3)
		real(dp) :: Sj(3), SiSj(3), Jij(3), eout
		real(dp), intent(inout) :: total_energy

		ion_ID = int(ion(0, i, j, k, l))
		! Over distinct bonds
		do nbdi = 1, no_of_nbd

			! For nbdi bond total connecting neighbours to the central ion [ID = ion(0, i, j, k, l)] is
			total_nbr = nn(nbdi, ion_ID, 0, 0)

			! no. of similar nbd for ith distinct bond (nbdi)
			do inbr = 1, total_nbr

				! nbr's cell position
				posx = nn(nbdi, ion_ID, inbr, 1)
				posy = nn(nbdi, ion_ID, inbr, 2)
				posz = nn(nbdi, ion_ID, inbr, 3)

				! nbr's ID in the cell
				cell = nn(nbdi, ion_ID, inbr, 4)
				ID_num = int(ion(4, posx, posy, posz, cell))
				if(ID_num.eq.0) then
					go to 4
				end if

				Sj(1:3) = ion(1:3, posx, posy, posz, cell)

				!Si.Sj
				SiSj = Si*Sj	!point-wise multiplication

				!Jij term
				Jij = j_exc(nbdi, central_ion_ID, ID_num, 1:3)

				!Jij.Si.Sj term
				eout = dot_product(SiSj, Jij)

				total_energy = total_energy + eout

4				continue
			end do
		end do

	end subroutine JSiSj
	
	!gmbSH	
	subroutine gmbSH(Si, g, mu, H, total_energy)

		implicit none

		real(dp) :: eout
		real(dp), intent(in) :: Si(3), g, mu, H(3)
		real(dp), intent(inout) :: total_energy

		! energy due to magnetic field
		eout = -(g*mu*dot_product(Si, H))
		
		total_energy = total_energy + eout
	
	end subroutine gmbSH

	!SINGLE ION ANISOTRPY (SIA)
	subroutine siaS2(Si, sia_val, total_energy)

		implicit none

		real(dp) :: Si2(3), eout
		real(dp), intent(in) :: Si(3), sia_val(3)
		real(dp), intent(inout) :: total_energy

		! Si**2
		Si2 = Si**2

		! energy due to single ion anisiatropy
		eout = dot_product(Si2, sia_val)			

		total_energy = total_energy + eout

	end subroutine siaS2

	end subroutine Hamiltonian
	
