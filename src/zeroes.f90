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
		real(dp) :: zero = real(0.0, dp)

		select case(zero_opt)

		case('sample')

			s_eng_avg = zero; e_eng_avg = zero
			s_eng2_avg = zero; e_eng2_avg = zero
			s_U_eng = zero; e_U_eng = zero
			s_cv = zero; e_cv = zero

			s_mag_avg = zero; e_mag_avg = zero
			s_mag2_avg = zero; e_mag2_avg = zero
			s_U_mag = zero; e_U_mag = zero
			s_chi = zero; e_chi = zero

		case('eng_mag')

			eng_avg = zero; eng2_avg = zero; eng4_avg = zero
			mag_avg = zero; mag2_avg = zero; mag4_avg = zero
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
