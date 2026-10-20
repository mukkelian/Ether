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

	subroutine Hamiltonian(trial, ih, S_trial, S_old, total_eng)

		use init, only: dp, ion, zeeman, SIA, XYZ
		implicit none

		integer, intent(in) :: ih
		integer :: ionID
		
		real(dp) :: S_diff(3), S2_diff(3)
		real(dp), intent(in) :: S_trial(3), S_old(3)	! central ion
		logical, intent(in) :: trial
		real(dp), intent(out) :: total_eng

		total_eng = 0.0_dp
		if(trial) then
			S_diff = S_trial - S_old
			S2_diff = S_trial**2 - S_old**2
		else
			S_diff = S_old
			S2_diff = S_old**2
		end if

		! Central ion's species ID
		ionID = int(ion(4, ih))

		call JSiSj(trial, ih, S_diff, ionID, total_eng)
		if(Zeeman) call gmbSH(S_diff, ionID, total_eng)
		if(SIA.and.XYZ) call siaS2(S2_diff, ionID, total_eng)

	end subroutine Hamiltonian
