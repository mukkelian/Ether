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

	subroutine mc(at_step, acceptance_count_)

		use init

		implicit none

		integer, intent(in) :: at_step
		integer, intent(out) :: acceptance_count_

		! Perform Monte Carlo
		call Monte_Carlo(acceptance_count_)

		after_equilibration_step: if((at_step.gt.tmcs_eq).and.(mod(real(at_step), real(to_cal)).eq.0))then

			total_calculations = total_calculations + 1
			eng = 0; net_mag = 0
			call get_tot_energy(eng)
			call get_tot_magnetisation(net_mag)

	        	call evaluate_observables('eng')
			call evaluate_observables('mag')

		end if after_equilibration_step

	end subroutine mc
