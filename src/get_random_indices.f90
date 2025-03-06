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

	subroutine get_random_indices(ith, jth, kth, lth)

	use init

	implicit none

	integer, intent(out) :: ith, jth, kth, lth

	real(dp) :: rn

        call get_random_num(fromx*real(1, dp), tox*real(1, dp), rn); ith = nint(rn)	
        if(ith.lt.fromx) then
        	ith = fromx
        elseif(ith.gt.tox)then
        	ith = tox
        end if

        call get_random_num(fromy*real(1, dp), toy*real(1, dp), rn); jth = nint(rn)
        if(jth.lt.fromy) then
        	jth = fromy
        elseif(jth.gt.toy)then
        	jth = toy
        end if

        call get_random_num(fromz*real(1, dp), toz*real(1, dp), rn); kth = nint(rn)
        if(kth.lt.fromz) then
        	kth = fromz
        elseif(kth.gt.toz)then
        	kth = toz
        end if
	
        call get_random_num(real(1, dp), lattice_per_unit_cell*real(1, dp), rn); lth = nint(rn)
        if(lth.lt.1) then
        	lth = 1
        elseif(lth.gt.(lattice_per_unit_cell))then
        	lth = lattice_per_unit_cell
        end if

	end subroutine get_random_indices
