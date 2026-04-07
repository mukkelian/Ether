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
        
	integer, parameter :: dp = 8

	integer :: tmcs, tmcs_eq, total_species_to_include, sc(3), repeat, repeati, &
		ovrr_steps, ovrr_MCS, optbeta, nspecies, total_ions_per_cell, &
		no_of_nbd, j_ID(2), similar_bonds, num_of_threads, total_ions, &
		lattice_per_unit_cell, nscan, &
		fromx, fromy, fromz, tox, toy, toz, total_calculations, itemp, &
		spin_file_ID, gss_ID, &
		seed_value, exchange_interval, total_info = 8, nbd_capacity, &
		spiral_capacity

        integer, allocatable, dimension(:) :: ions, tions, bblx, bbly, bblz, &
        	seed, included_species_ID, total_temperatures

        integer, allocatable, dimension(:,:) :: temp_range, ifs(:, :)

        integer, allocatable, dimension(:,:,:,:) :: nn

        real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)

        real(dp) :: temp, ht, lt, tint, dphi, to_cal, g_factor, J_para, &
        	lp(3), abc(3, 3), kb = 8.617333262d-5, &
        	h(3), anisotropy(3), mb, beta_critria, &
		nbd_finding_criteria, convert_to_rad = pi/180.0_dp, &
		beta, eng, eng_avg, eng2_avg, eng4_avg, mag_avg, mag2_avg, mag4_avg, &
		s_eng_avg, s_eng2_avg, e_eng2_avg, s_U_eng, s_cv, &
		s_mag_avg, s_mag2_avg, e_mag2_avg, s_U_mag, s_chi, &
		net_mag(3), err_U_mag, err_chi, err_mag_avg, err_U_eng, err_cv, &
		err_eng_avg, ovrr_para, SCabc(3, 3), total_energy, total_mag(3), &
		acceptance_counting, acceptance_count, to_angle = 180.0_dp/pi, &
		ss_dis, ss_direc(3), ss_latency, ss_proj(3), e_spiral_state, &
		spiral_state_avg, spiral_state2_avg, spiral_state4_avg, &
		s_spiral_state_avg, s_spiral_state2_avg, e_spiral_state2_avg, &
		s_U_spiral_state, s_spiral_state, err_U_spiral_state, err_spiral_state, & 
		err_spiral_state_avg, err_spiral_state_chi, s_spiral_state_chi

        real(dp), allocatable, dimension(:) :: x, y, z, s, nbd_dis, &
        	stgg, stgg_ion, &
        	e_mag_avg, e_U_mag, e_chi, e_eng_avg, e_U_eng, e_cv, &
		temp_T, s_mag_avg_T, s_chi_T, err_mag_avg_T, err_chi_T, &
		s_U_mag_T, err_U_mag_T, s_eng_avg_T, s_cv_T, err_eng_avg_T, &
		err_cv_T, s_U_eng_T, err_U_eng_T, &
		acceptance_ratio, temperature, s_spiral_state_avg_T, &
		s_spiral_state_chi_T, err_spiral_state_avg_T, &
		err_spiral_state_chi_T, s_U_spiral_state_T, err_U_spiral_state_T, &
		e_spiral_state_avg, e_U_spiral_state, e_spiral_state_chi

        real(dp), allocatable, dimension(:,:) :: sia_vec, mm_vector, ion
        
        real(dp), allocatable, dimension(:,:,:) :: mm_vector_avg_T

        real(dp), allocatable, dimension(:,:,:,:) :: j_exc

	character :: bc(3)
	character(len=2), allocatable, dimension(:) :: species_to_include, species
	character(len=5) :: model
	character(len=20), dimension(50) :: m_head, e_head, ss_head
	character(len=30) :: title, coordinate, lbl
	character(len=200), allocatable, dimension(:) :: input_data

	logical :: staggered, angle, Zeeman, SIA, para, ovrr, &
		PTalgo, temp_ex, beta_file, initiate_spin_files
		
!	For MPI's
	integer :: rank, size, ierr, total_observables, status(MPI_STATUS_SIZE), nprocs

	integer, allocatable, dimension(:) :: istart

	real(dp), allocatable, dimension(:) :: temp_assigned

	logical :: completed, Ising, XYZ, Checkerboard, ssp
	
        end module init
