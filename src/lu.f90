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

	!SMALL CAPS <----> LARGE CAPS
	subroutine lu(txtR, txtL, case_LU)

		implicit none

		character (len=*), intent(in) :: txtR
		character (len(txtR)), intent(out) :: txtL
        	character(len=54):: s
		character, intent(in) :: case_LU
		integer :: i, j

        	s = ' abcdefghijklmnopqrstuvwxyz ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		txtL = txtR
        	do i = 1, len(txtR)
        	        do j = 1, len(s)
        	                if( txtR(i:i).eq.s(j:j) )then
        	                        if( (j.le.27) .and. (case_LU.eq.'U') ) txtL(i:i) = s(j+27:j+27) ! into upper case
        	                        if( (j.gt.27) .and. (case_LU.eq.'L') ) txtL(i:i) = s(j-27:j-27) ! into lower case
        	                end if
        	        end do
        	end do
	end subroutine lu

