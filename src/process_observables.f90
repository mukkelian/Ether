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

	subroutine process_observables(observable_case)

		use init

		implicit none
		
		integer :: i, j, m
		
		character(len=*), intent(in) :: observable_case
		
		select case(observable_case)
		
		case('store')

			! ENERGY
			!temp_T(itemp) = temp; s_mag_avg_T(itemp) = s_mag_avg; s_chi_T(itemp) = s_chi
			!err_mag_avg_T(itemp) = err_mag_avg; err_chi_T(itemp) = err_chi; s_U_mag_T(itemp) = s_U_mag
			!err_U_mag_T(itemp) = err_U_mag

			local_obs(1 + li_obs) 	= temp
			local_obs(2 + li_obs) 	= s_mag_avg
			local_obs(3 + li_obs)	= s_chi
			local_obs(4 + li_obs) 	= err_mag_avg
			local_obs(5 + li_obs) 	= err_chi
			local_obs(6 + li_obs) 	= s_U_mag
			local_obs(7 + li_obs) 	= err_U_mag

		        ! MAGNETIC
!		        s_eng_avg_T(itemp) = s_eng_avg; s_cv_T(itemp) = s_cv; err_eng_avg_T(itemp) = err_eng_avg
!			err_cv_T(itemp) = err_cv; s_U_eng_T(itemp) = s_U_eng; err_U_eng_T(itemp) = err_U_eng

		        local_obs(8 + li_obs) 	= s_eng_avg
		        local_obs(9 + li_obs) 	= s_cv
		        local_obs(10 + li_obs) 	= err_eng_avg
			local_obs(11 + li_obs) 	= err_cv
			local_obs(12 + li_obs) 	= s_U_eng
			local_obs(13 + li_obs) 	= err_U_eng

			! ACCEPTANCE RATIOS
			local_obs(14 + li_obs) = acceptance_counting*100

			j = 0
			do i = 1, nspecies

				j = j + 1
				local_obs(total_observables + (3*(i-1) + 1) + li_obs) = mm_vector(j, 1)
				local_obs(total_observables + (3*(i-1) + 1) + 1 + li_obs) = mm_vector(j, 2)
				local_obs(total_observables + (3*(i-1) + 1) + 2 + li_obs) = mm_vector(j, 3)

			end do
			
			m = 0
			do i = 1, total_ions

				local_spn(m + 1 + li_spn) = ion(1, i)*s(int(ion(4, i)))
				local_spn(m + 2 + li_spn) = ion(2, i)*s(int(ion(4, i)))
				local_spn(m + 3 + li_spn) = ion(3, i)*s(int(ion(4, i)))
				local_spn(m + 4 + li_spn) = ion(4, i)
				m = m + 4

			end do

		case('write')

		lobs = 0	! local observable lenght
		lspn = 0	! local spin length
		do itemp = 1, nscan

			lobs = (itemp -1)*local_olen
			lspn = (itemp -1)*local_slen

			temp_T(itemp) 		= global_obs(1 + lobs)
			s_mag_avg_T(itemp) 	= global_obs(2 + lobs)
			s_chi_T(itemp) 		= global_obs(3 + lobs)
			err_mag_avg_T(itemp) 	= global_obs(4 + lobs)
			err_chi_T(itemp) 	= global_obs(5 + lobs)
			s_U_mag_T(itemp) 	= global_obs(6 + lobs)
			err_U_mag_T(itemp) 	= global_obs(7 + lobs)

			s_eng_avg_T(itemp) 	= global_obs(8 + lobs)
			s_cv_T(itemp) 		= global_obs(9 + lobs)
			err_eng_avg_T(itemp) 	= global_obs(10 + lobs)
			err_cv_T(itemp) 	= global_obs(11 + lobs)
			s_U_eng_T(itemp) 	= global_obs(12 + lobs)
			err_U_eng_T(itemp) 	= global_obs(13 + lobs)

			acceptance_ratio(itemp) = global_obs(14 + lobs)

			j = 0

			do i = 1, nspecies

			j = j + 1
			mm_vector_avg_T(j, 1, itemp) = global_obs(total_observables + (3*(i-1) + 1) + lobs)
			mm_vector_avg_T(j, 2, itemp) = global_obs(total_observables + (3*(i-1) + 1) + 1 + lobs)
			mm_vector_avg_T(j, 3, itemp) = global_obs(total_observables + (3*(i-1) + 1) + 2 + lobs)

			end do

			! Writing output files
			call write_output_files(itemp)
			! Writing Ground Spin States (GSS) file
			call write_gss(itemp)

		end do

		case default

			write(6, *) ''
			write(6, "(' ==> Found unknown case tag :',A8 )") observable_case
			write(6, *) "     in 'evaluate_observables' subroutine"
			write(6, *) '     STOPPING now'
			stop

		end select
		
	end subroutine process_observables
