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

	subroutine write_spins_at_K(K, T)

		use init

		implicit none

		integer, intent(in) :: T
		integer :: i, j
		
		real(dp), intent(in) :: K

		character(len=30) :: filename

                if(initiate_spin_files) then

			initiate_spin_files = .FALSE.
			!FORMATION OF FILES @ temp. K
			spin_file_ID = 2001 + T
			write(filename, 100) K
100		   	format( f9.4,'spK.xsf')
			open(file = adjustl(filename), unit = spin_file_ID)
			write(spin_file_ID, *) 'CRYSTAL'
			write(spin_file_ID, *) 'PRIMVEC'
			write(spin_file_ID, *) (abc(1,j)*(sc(1)), j= 1,3)
			write(spin_file_ID, *) (abc(2,j)*(sc(2)), j= 1,3)
			write(spin_file_ID, *) (abc(3,j)*(sc(3)), j= 1,3)
			write(spin_file_ID, *) 'CONVEC'
			write(spin_file_ID, *) (abc(1,j)*(sc(1)), j= 1,3)
			write(spin_file_ID, *) (abc(2,j)*(sc(2)), j= 1,3)
			write(spin_file_ID, *) (abc(3,j)*(sc(3)), j= 1,3)
			write(spin_file_ID, *) 'PRIMCOORD'
			write(spin_file_ID, *) total_ions, " 1"
			return

		end if

		do i = 1, total_ions
		if (ion(4, i).ne.0) then
			! SPINS
			write(spin_file_ID, 101) species(int(ion(4, i))), &
			ion(6:8, i), ion(1:3, i)
		end if
		end do

101	   	format(A5,1x, f13.7, 1x, f13.7, 1x, f13.7, 1x, f13.7, 1x, f13.7, 1x, f13.7)

		close(spin_file_ID)
	
	end subroutine write_spins_at_K
