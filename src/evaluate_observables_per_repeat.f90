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

	subroutine evaluate_observables_per_repeat

		use init

		implicit none
		
		integer :: i
		real(dp) :: mag_value, repeat_1

		repeat_1 = 1/(repeat*(repeat - real(1.0, dp)))
	
		!MAGNETIZATION PER REPEAT
		s_mag_avg = s_mag_avg/repeat
		s_U_mag = s_U_mag/repeat
		s_chi = s_chi/repeat

		! ERRORS
		err_U_mag = real(0, dp); err_chi = real(0, dp); err_mag_avg = real(0, dp)

		do i = 1, repeat

			err_U_mag = err_U_mag + (e_U_mag(i) - s_U_mag)**2
			err_chi = err_chi + (e_chi(i) - s_chi)**2
			err_mag_avg = err_mag_avg + (e_mag_avg(i) - s_mag_avg)**2

		end do

		err_U_mag = sqrt(repeat_1*err_U_mag)
		err_chi = sqrt(repeat_1*err_chi)
		err_mag_avg = sqrt(repeat_1*err_mag_avg)
	
		!ENERGY PER REPEAT
		s_eng_avg = s_eng_avg/repeat
		s_U_eng = s_U_eng/repeat
		s_cv = s_cv/repeat

		! ERRORS
		err_U_eng = real(0, dp); err_cv = real(0, dp); err_eng_avg = real(0, dp)

		do i = 1, repeat

			err_U_eng = err_U_eng + (e_U_eng(i) - s_U_eng)**2
			err_cv = err_cv + (e_cv(i) - s_cv)**2
			err_eng_avg = err_eng_avg + (e_eng_avg(i) - s_eng_avg)**2

		end do

		err_U_eng = sqrt(repeat_1*err_U_eng)
		err_cv = sqrt(repeat_1*err_cv)
		err_eng_avg = sqrt(repeat_1*err_eng_avg)
	
	end subroutine evaluate_observables_per_repeat
