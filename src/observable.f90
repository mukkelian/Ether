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

	subroutine observable(beta_value, at_step)

		use init

		implicit none

		integer, intent(in) :: at_step
		
		real(dp), intent(in) :: beta_value

		after_equilibration_step: if((at_step.gt.tmcs_eq).and.&
			(mod(real(at_step), real(to_cal)).eq.0))then

			total_calculations = total_calculations + 1
			eng = 0.0_dp; net_mag = 0.0_dp

			call get_tot_energy(eng)
			call get_tot_magnetisation(net_mag)

	        	call get_eng(beta_value, 'eng')
			call get_mag(beta_value, 'mag')
			if(ssp) call get_spiral_state(beta_value, 'spiral_state')

		end if after_equilibration_step

	end subroutine observable
