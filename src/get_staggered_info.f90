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

	allocate(stgg_ion(s_count), stg_IDs(nspecies))
	stg_IDs = 1

	inquire(file='staggered', exist=file_found)
	
	if(file_found.and.staggered) then

        	open(unit=2, file='staggered', status='old', action='read')
        	total_lines = 0
        	do
        		read(2, '(a)', end=10) lbl
        		total_lines = total_lines + 1
        	end do
10		continue
		close(2)

        	open(unit=2, file='staggered', status='old', action='read')

		allocate(stg_info(total_lines))

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
        end if

	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
	!$OMP DO SCHEDULE(DYNAMIC)
	do i = 1, total_ions
		stgg_ion(i) = stg_IDs(int(ion(4, i)))
	end do
	!$OMP END DO
	!$OMP END PARALLEL

	if (rank == 0) then
		write(6, *) ''
		write(6, *) '    :::::::::::::::::::::::::: &
        	                        STG list ::::::::::::::::::::::::'
		write(6, *) ''
        	if(staggered) write(6, *) "    (From 'staggered' file)"
        	if(.not.staggered) write(6, *) '    (DEFAULT values)'
        	write(6, *) ''
		do i = 1, nspecies
			write(6, "(9x,f7.1,4x,'-->',4x,A3)") &
				stg_IDs(i), species(i) 
			write(6, *) ''
		end do
	end if

	end subroutine get_staggered_info
