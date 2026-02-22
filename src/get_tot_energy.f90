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

	subroutine get_tot_energy(toten)

	use init
	use omp_lib

	implicit none

	integer :: i
	real(dp), intent(out) :: toten
	real(dp) :: partial_toten, Spin_vec(3)

	toten = 0.0
	call omp_set_nested(.true.)

	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, Spin_vec, partial_toten) REDUCTION(+:toten)
	!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(3)
	do i = 1, total_ions
			! Get spin vector from the ion array
			Spin_vec(1:3) = ion(1:3, i)

			! Calculate partial energy for this spin
			call Hamiltonian(.FALSE., i, Spin_vec, Spin_vec, partial_toten)

			! Accumulate the total energy
			toten = toten + partial_toten
		end do

			end do
		end do
	end do
	!$OMP END DO
	!$OMP END PARALLEL

	end subroutine get_tot_energy

