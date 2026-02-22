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

	subroutine Monte_Carlo(accept_count)

	use init
	use omp_lib

	implicit none

	integer :: i
	real(dp) :: S_vec_previous(5), S_vec_updated(5), &
		total_eng, eta
	integer, intent(out) :: accept_count

	accept_count = 0
	call omp_set_nested(.true.)

	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, &
	!$OMP& S_vec_previous, S_vec_updated, total_eng, eta)

	!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(3)
	do i = 1, total_ions

		! Copy current spin state
		S_vec_previous(1:5) = ion(1:5, i)

		! Update spin details
		call update_spin_details(S_vec_previous, S_vec_updated)

		! Calculate energy from Hamiltonian
		call Hamiltonian(.TRUE., i, S_vec_updated(1:3), &
				S_vec_previous(1:3), total_eng)

		! Get random number for Metropolis criterion
		call get_random_num(0d0, 1d0, eta)

		! Metropolis acceptance algorithm
		if (exp(-beta*total_eng) .gt. eta) then
		ion(1:5, i) = S_vec_updated(1:5)

		!$OMP ATOMIC
		accept_count = accept_count + 1
		!$OMP END ATOMIC

		end if

		end do

			end do
		end do
	end do

	!$OMP END DO
	!$OMP END PARALLEL

	end subroutine Monte_Carlo

