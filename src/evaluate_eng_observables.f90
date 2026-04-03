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

	subroutine evaluate_eng_observables(beta_value, observable_case)

		use init

		implicit none
		
		real(dp), intent(in) :: beta_value
		real(dp) :: U_eng, cv, eng_per_site
		
		character(len=*), intent(in) :: observable_case
		character (len=200) :: remark
		
		select case(observable_case)
		
		case('eng')
			
			eng_per_site = eng/(2*total_ions)
			eng_avg = eng_avg + eng_per_site
			eng2_avg = eng2_avg + eng_per_site**2
			eng4_avg = eng4_avg + eng_per_site**4

		case('avg_eng')

			eng_avg = eng_avg/total_calculations
			eng2_avg = eng2_avg/total_calculations
			eng4_avg = eng4_avg/total_calculations

			s_eng_avg = eng_avg + s_eng_avg
			e_eng_avg(repeati) = eng_avg

			U_eng = 1.0_dp - (1.0_dp/3.0_dp)*(eng4_avg/(eng2_avg**2))
			s_U_eng = U_eng + s_U_eng
			e_U_eng(repeati) = U_eng

			cv = (beta_value**2)*(eng2_avg - eng_avg**2)*total_ions
			s_cv = cv + s_cv
			e_cv(repeati) = cv

		case default

			remark = "Found unknown case tag '"//observable_case//&
			"' in evaluate_eng_observables"
			call terminate (remark)
		end select
		
	end subroutine evaluate_eng_observables
