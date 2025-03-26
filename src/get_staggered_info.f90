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

	subroutine get_staggered_info(s_count)
	
	use init
	
	implicit none
	
	integer, intent(in) :: s_count
	integer :: i
	logical :: file_found

	allocate(stgg_ion(s_count))
	
	inquire(file='staggered', exist=file_found)
	
	if(file_found) then
		if (rank == 0) write(6, *) ''
		if (rank == 0) write(6, *) '    :::::::::::::::::::::::::: &
        	                        STG list ::::::::::::::::::::::::'
        	open(unit=2, file='staggered', status='old', action='read')
        	do i = 1, s_count
        		read(2, *) stgg_ion(i)
        	end do
        	close(2)
		if (rank == 0) write(6, "(/,'      > ',8f7.1)") stgg_ion
		if (rank == 0) write(6, *) ''
	else
	if (rank == 0) then
		write(6, *) ''
		write(6, *) '	Logic for staggered magnetiaztion is .TRUE.'
		write(6, *) "	but 'staggered' file is not provided"
		write(6, *) '	     ~~~~~~~~~'
		write(6, *) ''
		write(6, *) "	file formate for 'staggered' is as follows"
		write(6, *) ''
		write(6, *) '	c1	for atom1'
		write(6, *) '	c2	for atom1'
		write(6, *) '	c3	for atom3'
		write(6, *) '	.		.'
		write(6, *) '	.		.'
		write(6, *) '	.		.'
		write(6, *) '	Note: Sequence of c[i] for ith atom should'
		write(6, *) '	       follows the sequnce as according to'
		write(6, *) '	       provided structure file'
		write(6, *) ''
		write(6, *) '	STOPPING now'
		write(6, *) ''
		stop
	end if
	end if

	end subroutine get_staggered_info
