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

	subroutine Hamiltonian(check_trail_spin, ih, S_trial, S_old, total_eng)
	
		use init

		implicit none

		integer, intent(in) :: ih
		integer :: ionID, nbdi, inbr, cell
		
		real(dp) :: sia_value(3), S_diff(3), S2_diff(3)
		real(dp), intent(in) :: S_trial(3), S_old(3)	! central ion
		logical, intent(in) :: check_trail_spin
		real(dp), intent(out) :: total_eng

		total_eng = 0
		if(check_trail_spin) then
			S_diff = S_trial - S_old
			S2_diff = S_trial**2 - S_old**2
		else
			S_diff = S_old
			S2_diff = S_old**2
		end if

		! Central ion's species ID
		ionID = int(ion(4, ih))

		call JSiSj(ih, S_diff, ionID, total_eng)

		if(Zeeman) call gmbSH(S_diff, g_factor, mb, H, total_eng)

		if(SIA.and.XYZ) then
			sia_value(1:3) = sia_vec(1:3, ionID)
			call siaS2(S2_diff, sia_value, total_eng)
		end if
	
	contains

	! JSiSj	
	subroutine JSiSj(i, Si, central_ion_ID, total_energy)
	
		implicit none

		integer, intent(in) :: i, central_ion_ID
		integer :: total_nbr, ID_num, ion_ID

		real(dp), intent(in) :: Si(3)
		real(dp) :: Sj(3), SiSj(3), Jij(3), eout
		real(dp), intent(inout) :: total_energy

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
	subroutine siaS2(Si2, sia_val, total_energy)

		implicit none

		real(dp) :: eout
		real(dp), intent(in) :: Si2(3), sia_val(3)
		real(dp), intent(inout) :: total_energy

		! energy due to single ion anisiatropy
		eout = dot_product(Si2, sia_val)

		total_energy = total_energy + eout

	end subroutine siaS2

	end subroutine Hamiltonian
	
