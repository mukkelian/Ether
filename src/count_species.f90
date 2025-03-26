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

	subroutine count_species(line, n)

	implicit none

	integer :: i, j, k
	integer, intent(out) :: n
	character(len=*), intent(in) :: line

		j = 0; k = 0

        i = 1
100     continue
        check: if (line(i:i).ne.' ') then
			k = i+1
102     	continue
        	if(line(k:k) .ne.' ')then
        		k = k+1
        	    goto 102
        	else
        	    j = j + 1
        	end if
        	i = k + 1
		else
				i = i + 1
		end if check

		if (i.lt.(len_trim(line)+1)) goto 100
		n = j	!no. of species

	end subroutine count_species