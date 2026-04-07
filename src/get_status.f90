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

	subroutine get_status(rankID, step_mcs, step_repeat, T)
	
	use init, only : dp, tmcs, repeat
	
	implicit none

	integer, intent(in) :: rankID, step_mcs, step_repeat
	real(dp), intent(in) :: T
	
	character(len=20) :: step1, step2, step3, step4

	write(step1, '(i0)') step_mcs
	write(step2, '(i0)') tmcs
	write(step3, '(i0)') step_repeat
	write(step4, '(i0)') repeat

	write(6, 1) rankID, T, trim(step3), trim(step4), &
            trim(step1), trim(step2)
	
1	format('STATUS >> Rank ID:',1x, i3, 4x,'Temperature:',1x, g11.4,1x,&
	'Repeat:',1x, A,'/', A, 1x,&
	'MSCS:',1x, A,'/',A)

	end subroutine get_status
