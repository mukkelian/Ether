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

	subroutine overrelaxation

	use init

	implicit none

	integer :: i, j, k, l, m

	real(8) :: A_ovr(3), Si(3), rn_ovr

        overrelaxation_method : do m = 1, overr_steps
		! Duplicate Names Between Modules are not Allowed (DNBM)
		! it causes variables to exchange between the modules/threads
	        call get_random_num(fromx*1d0, tox*1d0, rn_ovr); i = int(rn_ovr)	
	        if(i.lt.fromx) then
	        	i = fromx
	        elseif(i.gt.tox)then
	        	i = tox
	        end if

	        call get_random_num(fromy*1d0, toy*1d0, rn_ovr); j = int(rn_ovr)	!DNBM
	        if(j.lt.fromy) then
	        	j = fromy
	        elseif(j.gt.toy)then
	        	j = toy
	        end if

	        call get_random_num(fromz*1d0, toz*1d0, rn_ovr); k = int(rn_ovr)	!DNBM
	        if(k.lt.fromz) then
	        	k = fromz
	        elseif(k.gt.toz)then
	        	k = toz
	        end if
	
	        call get_random_num(1d0, lattice_per_unit_cell*1d0, rn_ovr); l = int(rn_ovr)	!DNBM
	        if(l.lt.1) then
	        	l = 1
	        elseif(l.gt.(lattice_per_unit_cell))then
	        	l = lattice_per_unit_cell
	        end if

		A_ovr = 0

		call get_ovrr_vec(i, j, k, l, A_ovr)
	

		Si(1:3) = ion(1:3, i, j, k, l)

	        Si(1:3) = (2* &
	        (dot_product(A_ovr(1:3), Si(1:3))/dot_product(A_ovr(1:3), A_ovr(1:3)) )* &
	        A_ovr(1:3)) - Si(1:3)

		ion(1:3, i, j, k, l) = Si(1:3)/sqrt(dot_product(Si(1:3), Si(1:3)))

	end do overrelaxation_method

	end subroutine overrelaxation
