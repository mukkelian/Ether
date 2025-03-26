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

	subroutine get_sia_values

	use init

	implicit none

	integer :: i, j, SIA_ID
	real(dp) :: sia_vectors(3)
	character(len=2) :: atom1, ab

	logical :: file_found

	if(sia.and.XYZ) then

	if(allocated(sia_vec)) deallocate(sia_vec)
	allocate(sia_vec(1:3, nspecies))
	
	inquire(file='single_ion_anisotropy', exist=file_found)
	if(file_found) then
		open(10011, file='single_ion_anisotropy', status='old', action='read')
		if (rank == 0) write(6, '(3X,":::::::::: SIA list (meV) ::::::::::")')
		if (rank == 0) write(6, *) ''
		do i = 1, nspecies
			read(10011, *) atom1, sia_vectors(1:3)	! Single Ion Anisotropy (SIA) vector in unit of meV
			sia_vec(1:3, i) = sia_vectors(1:3)
			ab = atom1
			call lu(ab(1:1), ab(1:1), 'U' )
			call lu(ab(2:2), ab(2:2), 'L' )
			atom1 = ab

			if (rank == 0) write(6, "(4X,A4,3f10.3)") atom1, sia_vec(1:3, i)

			do j = 1, nspecies
				if (atom1.eq.species(j)) then
					SIA_ID = j
				end if
			end do
			sia_vec(1:3, SIA_ID) = sia_vec(1:3, SIA_ID)*s(SIA_ID)**2
		end do
		write(6, *) ''
		close(10011)
	else
		if (rank == 0) then
			write(6, *) "	 SIA is .TRUE. but file 'single_ion_anisotropy'"
			write(6, *) "	 is not found so, stopping now"
			write(6, *) ''
			write(6, *) "	 'single_ion_anisotropy' file formate are as follows"
			write(6, *) "	  ~~~~~~~~~~~~~~~~~~~~~"
			write(6, *) "	 species1_name	1_SIAx	1_SIAy	1_SIAz"
			write(6, *) "	 species2_name	2_SIAx	2_SIAy	2_SIAz"
			write(6, *) "	 ."
			write(6, *) "	 ."
			write(6, *) "	 ."
			write(6, *) ''
		end if
		stop
	end if

	end if

	end subroutine get_sia_values
