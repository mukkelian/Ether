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

        module init

        use mpi

        implicit none

	integer, public, save :: tmcs, tmcs_eq, n_speci_incl, sc(3), sample, samplei, &
		overr_steps, ovrr_start_interval, optbeta, nspecies, total_ions_per_cell, &
		no_of_nbd, j_ID(2), similar_bonds, num_of_threads, total_lattice_sites, &
		lattice_per_unit_cell, nscan, nbd_cell_x, nbd_cell_y, nbd_cell_z, &
		total_ions, fromx, fromy, fromz, tox, toy, toz, total_calculations, itemp, &
		spin_file_ID, gss_ID, acceptance_counting, acceptance_count
	
	real(8), public, parameter :: pi = 4*atan(1d0)

        real(8), public, save :: temp, ht, lt, tint, dphi, to_cal, g_factor, para_value, &
        	exchange_interval, lp(3), abc(3, 3), sia(1:3), kb = 8.617333262d-5, &
        	h(3), anisotropy(3), j_value, mb, beta_critria, &
		nbd_finding_criteria = 0.0001d0, convert_to_rad = pi/180.0d0, beta, &
		eng, eng_avg, eng2_avg, eng4_avg, mag_avg, mag2_avg, mag4_avg, &
		s_eng_avg, s_eng2_avg, e_eng2_avg, s_U_eng, s_cv, &
		s_mag_avg, s_mag2_avg, e_mag2_avg, s_U_mag, s_chi, &
		net_mag(3), err_U_mag, err_chi, err_mag_avg, err_U_eng, err_cv, &
		err_eng_avg, overr_para

        integer, allocatable, public, save :: ionn(:), nn(:,:,:,:), bblx(:), bbly(:), bblz(:)

        real(8), allocatable, public, save :: x(:), y(:), z(:), s(:), nbd_dis(:), &
        	j_exc(:,:,:,:), sia_vec(:, :), stgg(:), stgg_ion(:), ion(:, :, :, :, :), &
        	e_mag_avg(:), e_U_mag(:), e_chi(:), e_eng_avg(:), e_U_eng(:), e_cv(:), &
        	mm_vector(:, :), &
		temp_T(:), s_mag_avg_T(:), s_chi_T(:), err_mag_avg_T(:), err_chi_T(:), &
		s_U_mag_T(:), err_U_mag_T(:), s_eng_avg_T(:), s_cv_T(:), err_eng_avg_T(:), &
		err_cv_T(:), s_U_eng_T(:), err_U_eng_T(:), mm_vector_avg_T(:, :, :), &
		acceptance_ratio(:)

	character, public, save :: bc(3)
	character(len=2), allocatable, public, save :: specicies_to_include(:), species(:)
	character(len=5), public, save :: model
	character(len=20), dimension(50), public, save :: m_head, e_head
	character(len=30), public, save :: title, coordinate, filename, lbl

	logical, public, save :: staggered, angle, Zeeman, single_ion_anisotropy, para, ovrr, &
		EXalgo, temp_ex, beta_file, initiate_spin_files
		
!	For MPI's
	integer :: rank, size, ierr, interval, left, tag, local_olen, local_slen, &
		total_observables, status(MPI_STATUS_SIZE), request_obs, request_spn, &
		status_obs, status_spn, li_obs, li_spn, lobs, lspn

	integer, allocatable :: addmoreitr(:), num_iterations(:), &
		istart(:), iend(:), tlobs(:), tlspn(:), si_obs(:), &
		si_spn(:)

	real(8), allocatable, public, save :: local_obs(:), local_spn(:), &
		global_obs(:), global_spn(:), temp_assigned(:)

	logical, public, save :: completed, Ising, XYZ
	
        end module init
