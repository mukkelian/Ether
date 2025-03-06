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

	subroutine j_values
	
	use init

	implicit none
	
	integer :: i, j, k, ii, nbd_dis_max

	character(len=3) :: atom(2), atom1, atom2, ab, out
	logical :: file_present
		
	inquire(file='j_exchange', exist=file_present)
	if(.not.file_present) then
	if (rank == 0) then
		write(6, *) "==> Required file 'j_exchange' is not present"
		write(6, *) "	 STOPPING now"
		write(6, *) ""
	end if
		stop
	end if

	open(2, file='j_exchange', status='old', action='read')
	read(2, *) no_of_nbd					! no. of diff. bond length w.r.t. any central ion
	if(allocated(nbd_dis)) deallocate(nbd_dis)
	allocate(nbd_dis(no_of_nbd))				! nbd distances
	if(allocated(j_exc)) deallocate(j_exc)
	allocate(j_exc(no_of_nbd, nspecies, nspecies, 1:3))    	! j_exc(:,:,xx-yy-zz)
	j_exc = real(0.0, dp); j_ID = 0

	if (rank == 0) then
		write(6, *) '==> J-values are recieved'
		write(6, *) ''
		write(6, *) '  ::::::::: J values (meV) ::::::::::::'
		write(6, *) ''
		write(6, *) '  ion1 <---> ion2     J    distance (A)'
	end if

	do i = 1, no_of_nbd
		read(2, *) nbd_dis(i), similar_bonds
		do ii = 1, similar_bonds

			read(2, *) atom1, atom2, j_value, anisotropy(1:3)
			ab = atom1
			call lu(ab(1:1), out, 'U' )
			atom1(1:1) = out
			call lu(ab(2:2), out, 'L' )
			atom1(2:2) = out
		
			ab = atom2	
			call lu(ab(1:1), out, 'U' )
			atom2(1:1) = out
			call lu(ab(2:2), out, 'L' )
			atom2(2:2) = out
			atom = (/ atom1, atom2 /)

			do j = 1, 2
				do k = 1, nspecies
					if (atom(j).eq.species(k)) then
						j_ID(j) = k
					end if
				end do
			end do
			if (rank == 0) write(6, 501) species(j_ID(1)), species(j_ID(2)), j_value, &
			nbd_dis(i)
501			format(5X,A2,X,'<--->',X,A2,2X,f8.3,f9.5)
			j_exc(i, j_ID(1), j_ID(2), 1:3) = j_value*anisotropy*s(j_ID(1))*s(j_ID(2))
			j_exc(i, j_ID(2), j_ID(1), 1:3) = j_exc(i, j_ID(1), j_ID(2), 1:3)

		end do
	end do

	!closing of j_exchange file
	close(2)

	if (rank == 0) then
	
		write(6, *) ''
		write(6, *) '    Note: Kindly check the above listed J values'
		write(6, *) '    ````` from j_exchange file before proceeding further!'
		write(6, *) ''

	end if

	end subroutine j_values
