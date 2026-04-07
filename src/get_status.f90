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
	
	integer, dimension(8) :: time_stamp	
	integer, intent(in) :: rankID, step_mcs, step_repeat
	real(dp), intent(in) :: T
	
	write(6, 1) rankID, T, step_repeat, repeat,step_mcs,tmcs, &
	time_stamp(3),time_stamp(2),time_stamp(1), time_stamp(5:7)

1	format('Rank ID:',1x, i3, 4x,'STATUS --> Temperature:',1x, g11.4,1x,&
	'Repeat:',1x, i2, 1x,'/',1x, i2, 1x,&
	'MSCS:',1x, i7,1x,'/',1x,i7,/,10x,&
	'on date'1x,i2,'-',i2,'-',i4,4x,&
	'at time',1x,i2,1x,'hrs.',1x,i2,1x,'min.',1x,i2,1x,'sec.',//)

	end subroutine get_status
