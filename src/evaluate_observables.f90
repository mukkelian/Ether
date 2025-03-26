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

	subroutine evaluate_observables(observable_case)

		use init

		implicit none
		
		integer :: i, j, k, l, m
		real(dp) :: mag_value, U_mag, chi, U_eng, cv, repeat_1, mag_per_site, &
			eng_per_site
		
		character(len=*), intent(in) :: observable_case
		
		select case(observable_case)
		
		case('eng')
			
			eng_per_site = eng/(2*total_lattice_sites)
			eng_avg = eng_avg + eng_per_site
			eng2_avg = eng2_avg + eng_per_site**2
			eng4_avg = eng4_avg + eng_per_site**4

		case('mag')

			mag_value =  sqrt(dot_product(net_mag, net_mag))
			mag_per_site = mag_value/total_lattice_sites
			mag_avg = mag_avg + mag_per_site
			mag2_avg = mag2_avg + mag_per_site**2
			mag4_avg = mag4_avg + mag_per_site**4

		case('avg_mag')

			mag_avg = mag_avg/total_calculations
			mag2_avg = mag2_avg/total_calculations
			mag4_avg = mag4_avg/total_calculations
			s_mag_avg = mag_avg + s_mag_avg
			e_mag_avg(repeati) = mag_avg

			U_mag = real(1, dp) - ((real(1, dp)/real(3, dp))*(mag4_avg/(mag2_avg**2)))
			s_U_mag = U_mag + s_U_mag
			e_U_mag(repeati) = U_mag
			chi = beta*(mag2_avg - mag_avg**2)*total_lattice_sites
			s_chi = chi + s_chi
			e_chi(repeati) = chi

		case('avg_eng')

			eng_avg = eng_avg/total_calculations
			eng2_avg = eng2_avg/total_calculations
			eng4_avg = eng4_avg/total_calculations

			s_eng_avg = eng_avg + s_eng_avg
			e_eng_avg(repeati) = eng_avg

			U_eng = 1 - ((real(1, dp)/real(3, dp))*(eng4_avg/(eng2_avg**2)))
			s_U_eng = U_eng + s_U_eng
			e_U_eng(repeati) = U_eng

			cv = (beta**2)*(eng2_avg - eng_avg**2)*total_lattice_sites
			s_cv = cv + s_cv
			e_cv(repeati) = cv

		case('observables_per_repeat')

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
			local_obs(14 + li_obs) = acceptance_counting/&
				(tmcs*total_lattice_sites*repeat*real(1, dp))*100

			j = 0
			do i = 1, nspecies

				j = j + 1
				local_obs(total_observables + (3*(i-1) + 1) + li_obs) = mm_vector(j, 1)
				local_obs(total_observables + (3*(i-1) + 1) + 1 + li_obs) = mm_vector(j, 2)
				local_obs(total_observables + (3*(i-1) + 1) + 2 + li_obs) = mm_vector(j, 3)

			end do
			
			m = 0
			do i = fromx, tox
				do j = fromy, toy
					do k = fromz, toz
						do l = 1, lattice_per_unit_cell

				local_spn(m + 1 + li_spn) = ion(1, i, j, k, l)*s(int(ion(4, i, j, k, l)))
				local_spn(m + 2 + li_spn) = ion(2, i, j, k, l)*s(int(ion(4, i, j, k, l)))
				local_spn(m + 3 + li_spn) = ion(3, i, j, k, l)*s(int(ion(4, i, j, k, l)))
				local_spn(m + 4 + li_spn) = ion(4, i, j, k, l)
				m = m + 4

						end do
					end do
				end do
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
		
	end subroutine evaluate_observables
