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

	subroutine zeroes(zero_opt)
	
		use init

		implicit none
		
		character(len=*), intent(in) :: zero_opt

		select case(zero_opt)

		case('sample')

			s_eng_avg = 0.0; e_eng_avg = 0.0
			s_eng2_avg = 0.0; e_eng2_avg = 0.0
			s_U_eng = 0.0; e_U_eng = 0.0
			s_cv = 0.0; e_cv = 0.0

			s_mag_avg = 0.0; e_mag_avg = 0.0
			s_mag2_avg = 0.0; e_mag2_avg = 0.0
			s_U_mag = 0.0; e_U_mag = 0.0
			s_chi = 0.0; e_chi = 0.0

		case('eng_mag')

			eng_avg = 0.0; eng2_avg = 0.0; eng4_avg = 0.0
			mag_avg = 0.0; mag2_avg = 0.0; mag4_avg = 0.0
			total_calculations = 0

		case default

		if (rank == 0) then
			write(6, *) ''
			write(6, "(' ==> Found unknown case tag :',A8 )") zero_opt
			write(6, *) "in 'zeroes' subroutine"
			write(6, *) 'STOPPING now'
			stop
		end if

		end select


	end subroutine zeroes
