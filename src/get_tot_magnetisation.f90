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

	subroutine get_tot_magnetisation(magnetisation)

	use init
	use omp_lib

	implicit none

	integer :: i, stg
	real(dp), intent(out) :: magnetisation(3)
	real(dp) :: temp_magnetisation(3)

	magnetisation = real(0, dp)
	call omp_set_nested(.true.)

	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, stg, temp_magnetisation) &
	!$OMP& REDUCTION(+:magnetisation)  ! Reduction on all elements of the magnetisation array

	temp_magnetisation = real(0, dp)

	!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(3)
	do i = 1, total_ions

		! Get the stg value depending on the 'staggered' flag
		if (staggered .eqv. .TRUE.) then
			stg = stgg_ion(l)
		else
			stg = 1
		end if

		temp_magnetisation(1:3) = temp_magnetisation(1:3) + &
			stg*ion(1:3, i)*s(int(ion(4, i)))

	end do
	!$OMP END DO

	! Accumulating net magnetisation
	magnetisation(1:3) = magnetisation(1:3) + temp_magnetisation(1:3)

	!$OMP END PARALLEL  ! End of the parallel region

	end subroutine get_tot_magnetisation

