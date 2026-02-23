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
	use omp_lib

	implicit none
	
	integer, intent(in) :: s_count
	integer :: i, k, total_lines
	
	real(dp), allocatable :: stg_info(:), stg_IDs(:)
	logical :: file_found
	
	character(len=3) :: atom, ab, out

	allocate(stgg_ion(s_count))
	
	inquire(file='staggered', exist=file_found)
	
	if(file_found) then
		if (rank == 0) write(6, *) ''
		if (rank == 0) write(6, *) '    :::::::::::::::::::::::::: &
        	                        STG list ::::::::::::::::::::::::'
        	open(unit=2, file='staggered', status='old', action='read')
        	total_lines = 0
        	do
        		read(2, '(a)', end=10) lbl
        		total_lines = total_lines + 1
        	end do
10		continue
		close(2)

        	open(unit=2, file='staggered', status='old', action='read')

		allocate(stg_info(total_lines), stg_IDs(nspecies))
		stg_IDs = 1
        	do i = 1, total_lines
        		read(2, *) stg_info(i), atom

			ab = trim(adjustl(atom))
			call lu(ab(1:1), out, 'U' )
			atom(1:1) = out
			call lu(ab(2:2), out, 'L' )
			atom(2:2) = out

			do k = 1, nspecies
				if (atom.eq.species(k)) then
					stg_IDs(k) = stg_info(i)
				end if
			end do
			
        	end do
        	close(2)

		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
		!$OMP DO SCHEDULE(DYNAMIC)
		do i = 1, total_ions
			stgg_ion(i) = stg_IDs(int(ion(4, i)))
		end do
		!$OMP END DO
		!$OMP END PARALLEL

		if (rank == 0) then
			do i = 1, nspecies
				write(6, "(9x,8f7.1,4x,'-->'4x,A3)") &
					stg_IDs(i), species(i) 
				write(6, *) ''
			end do
		end if
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
