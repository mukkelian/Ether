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

	subroutine evaluate_mag_observables(beta_value, observable_case)

		use init

		implicit none

		real(dp), intent(in) :: beta_value
		real(dp) :: mag_value, U_mag, chi, mag_per_site
		
		character(len=*), intent(in) :: observable_case
		character(len=200) :: remark
		
		select case(observable_case)

		case('mag')

			mag_value =  sqrt(dot_product(net_mag, net_mag))
			mag_per_site = mag_value/total_ions
			mag_avg = mag_avg + mag_per_site
			mag2_avg = mag2_avg + mag_per_site**2
			mag4_avg = mag4_avg + mag_per_site**4

		case('avg_mag')

			mag_avg = mag_avg/total_calculations
			mag2_avg = mag2_avg/total_calculations
			mag4_avg = mag4_avg/total_calculations
			s_mag_avg = mag_avg + s_mag_avg
			e_mag_avg(repeati) = mag_avg

			U_mag = 1.0_dp - (1.0_dp/3.0_dp)*(mag4_avg/(mag2_avg**2))
			s_U_mag = U_mag + s_U_mag
			e_U_mag(repeati) = U_mag
			chi = beta_value*(mag2_avg - mag_avg**2)*total_ions
			s_chi = chi + s_chi
			e_chi(repeati) = chi

		case default

			remark = "Found unknown case tag '"//observable_case//&
			"' in evaluate_mag_observables"
			call terminate(remark)

		end select
		
	end subroutine evaluate_mag_observables
