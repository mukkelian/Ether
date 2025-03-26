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

	subroutine allocate_observables

		use init

		implicit none

		allocate(e_mag_avg(repeat), e_U_mag(repeat), e_chi(repeat), &
		e_eng_avg(repeat), e_U_eng(repeat), e_cv(repeat), acceptance_ratio(nscan), &
		mm_vector(nspecies, 1:3), temp_T(nscan), s_mag_avg_T(nscan), s_chi_T(nscan), &
		err_mag_avg_T(nscan), err_chi_T(nscan), s_U_mag_T(nscan), err_U_mag_T(nscan), &
		s_eng_avg_T(nscan), s_cv_T(nscan), err_eng_avg_T(nscan), err_cv_T(nscan), &
		s_U_eng_T(nscan), err_U_eng_T(nscan), mm_vector_avg_T(nspecies, 1:3, nscan))

	end subroutine allocate_observables
