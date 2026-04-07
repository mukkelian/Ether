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

	subroutine parallel_tempering(ion, energy, &
		beta, itemp_num, swap_count, comm)

	use init, only: dp, total_ions, total_info
	use mpi

	implicit none

	integer, intent(in) :: swap_count, comm, itemp_num
	integer :: partner, ierr, count, &
		accept_int, local_rank, local_size

	real(dp), intent(inout) :: ion(0:total_info, total_ions)
	real(dp), intent(in) :: energy, beta
	real(dp) :: energy_p, beta_p, ion_p(0:total_info, total_ions), &
		delta, r

	logical :: accept

	call MPI_Comm_rank(comm, local_rank, ierr)
	call MPI_Comm_size(comm, local_size, ierr)

	accept = .false.

	! Skip dummy temperature
	if (itemp_num == 0) return

	count = 9*total_ions

	! Partner selection (odd-even)
	if (mod(swap_count,2) == 0) then
		if (mod(local_rank,2) == 0) then
			partner = local_rank + 1
		else
			partner = local_rank - 1
		end if
	else
		if (mod(local_rank,2) == 1) then
			partner = local_rank + 1
		else
			partner = local_rank - 1
		end if
	end if

	! Safety
	if (partner < 0 .or. partner >= local_size) return

	! Exchange energy and beta (ALL ranks)
	call MPI_Sendrecv(energy, 1, MPI_DOUBLE_PRECISION, partner, 1, &
		energy_p, 1, MPI_DOUBLE_PRECISION, partner, 1, &
		comm, MPI_STATUS_IGNORE, ierr)

	call MPI_Sendrecv(beta, 1, MPI_DOUBLE_PRECISION, partner, 2, &
		beta_p, 1, MPI_DOUBLE_PRECISION, partner, 2, &
		comm, MPI_STATUS_IGNORE, ierr)

	! Only lower rank computes
	if (local_rank < partner) then

		delta = -(beta_p - beta) * (energy - energy_p)

		if (delta > 0.0_dp) then
			accept = .true.
		else
			call get_random_num(0.0_dp, 1.0_dp, r)
			accept = (r < exp(delta))
		end if
	else
		accept = .false.
	end if

	! Share decision with partner
	accept_int = merge(1, 0, accept)

	call MPI_Sendrecv(accept_int, 1, MPI_INTEGER, partner, 3, &
		accept_int, 1, MPI_INTEGER, partner, 3, &
		comm, MPI_STATUS_IGNORE, ierr)

	accept = (accept_int == 1)

	! ALWAYS exchange spins
	call MPI_Sendrecv(ion, count, MPI_DOUBLE_PRECISION, partner, 4, &
		ion_p, count, MPI_DOUBLE_PRECISION, partner, 4, &
		comm, MPI_STATUS_IGNORE, ierr)

	! Apply swap
	if (accept) then
		ion = ion_p
	end if

	end subroutine parallel_tempering
