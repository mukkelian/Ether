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

	subroutine get_temperature()
	
	use init, only: temperature, ht, lt, tint, AM, nscan, root
	
	implicit none

	integer :: i
	character(len=2) :: lbl
	logical :: file_present

	file_present = .FALSE.

	if(am) inquire(file='grid.out', exist=file_present)
	if(am.and.file_present) then
		open(unit=11, file='grid.out', status='old', action='read')
		read(11, *) lbl, nscan
		if(allocated(temperature)) deallocate(temperature)
		allocate(temperature(nscan))
		read(11, *) lbl, temperature(1:nscan)
		ht = temperature(1)
		lt = temperature(nscan)
		if(root) then
		write(6, *) "==> New temperature grids from 'grid.out' [generated from &
			Adapted Mesh Tool] has been read"
		write(6, *) '	New temperature grids are:'
		write(6, '(8x,g11.4,1x,g11.4,1x,g11.4,1x,g11.4, 1x,g11.4)') temperature(1:nscan)
		end if
	else
		if(root) then
		write(6, *) "==> Adaptive mesh file 'grid.out' is absent"
		write(6, *) "	Continuing with equispaced temperature points as per descriptions present in 'in.ether'"
		end if
	end if

	if(.not.file_present) then
		! Total temperature count
		nscan = nint((ht - lt) / (tint)) + 1
		if(allocated(temperature)) deallocate(temperature)
		allocate(temperature(nscan))
		do i = 1, nscan
			temperature(i) =  ht - tint * (i - 1)
		end do
	end if

	end subroutine get_temperature
