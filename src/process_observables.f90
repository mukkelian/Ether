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

	subroutine process_observables(observable_case, T)

		use init

		implicit none
		
		integer :: i, j, k, n
		integer, intent(in) :: T
		
		real(dp), allocatable :: observables(:)
		
		character(len=*), intent(in) :: observable_case
		character (len=200) :: remark
		
		n = total_observables + 3*nspecies
		allocate(observables(n)) 
		select case(observable_case)
		
		case('store')

			observables(1) = temp
			observables(2) = s_mag_avg
			observables(3) = s_chi
			observables(4) = err_mag_avg
			observables(5) = err_chi
			observables(6) = s_U_mag
			observables(7) = err_U_mag
			observables(8) = s_eng_avg
			observables(9) = s_cv
			observables(10) = err_eng_avg
			observables(11) = err_cv
			observables(12) = s_U_eng
			observables(13) = err_U_eng
			observables(14) = acceptance_counting*100

			j = 0
			do i = 1, nspecies
				j = j + 1
				observables(total_observables + (3*(i-1) + 1)) = mm_vector(j, 1)
				observables(total_observables + (3*(i-1) + 1) + 1) = mm_vector(j, 2)
				observables(total_observables + (3*(i-1) + 1) + 2) = mm_vector(j, 3)
			end do

				
			! Writing oservables into ETHER.obs
        		call rw_file('w', 'ETHER.obs', (/1, n/), (/1, 1/), T, observables)

		case('write')
		
		! Loop over all temperature ID
		do k = 1, nscan

			! Reading oservables from ETHER.obs
			call rw_file('r', 'ETHER.obs', (/1, n/), (/1, 1/), k, observables)

			temp_T(k) = observables(1)
			s_mag_avg_T(k) 	= observables(2)
			s_chi_T(k) = observables(3)
			err_mag_avg_T(k) = observables(4)
			err_chi_T(k) = observables(5)
			s_U_mag_T(k) = observables(6)
			err_U_mag_T(k) = observables(7)
			s_eng_avg_T(k) 	= observables(8)
			s_cv_T(k) = observables(9)
			err_eng_avg_T(k) = observables(10)
			err_cv_T(k) = observables(11)
			s_U_eng_T(k) = observables(12)
			err_U_eng_T(k) = observables(13)
			acceptance_ratio(k) = observables(14)

			j = 0

			do i = 1, nspecies

			j = j + 1
			mm_vector_avg_T(j, 1, k) = observables(total_observables + (3*(i-1) + 1))
			mm_vector_avg_T(j, 2, k) = observables(total_observables + (3*(i-1) + 1) + 1)
			mm_vector_avg_T(j, 3, k) = observables(total_observables + (3*(i-1) + 1) + 2)

			end do

			! Writing output files
			call write_output_files(k)
			! Writing Ground Spin States (GSS) file
			call write_gss(k)

		end do

		case default
			remark = "Found unknown case tag '"//observable_case//&
			"' in process_observables"
			call terminate(remark)
		end select
		
		deallocate(observables)
		
	end subroutine process_observables
