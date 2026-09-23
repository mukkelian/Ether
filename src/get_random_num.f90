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

	subroutine get_random_num(from, to, rn)

		use init, only : dp, random_algo
		use mt19937

		implicit none

		real(dp), intent(in) :: from, to
		real(dp), intent(out) :: rn

		select case(trim(adjustl(random_algo)))
		
		case('xoshiro256')
			call random_number(rn)
		case('mt19937')
			rn = mt_rand()
		case default
			call terminate("Found unknown/missing information: "//trim(adjustl(random_algo)))
		end select
		rn = (to-from)*rn + from

	end subroutine

	subroutine set_seed(offset)

		use init, only : seed, random_algo
		use mt19937

		implicit none

		integer, intent(in) :: offset

		select case(trim(adjustl(random_algo)))

		case('xoshiro256')
			call random_seed(put = seed + offset)
		case('mt19937')
			call mt_init(seed(1) + offset)
		case default
			call terminate("Found unknown/missing information: "//trim(adjustl(random_algo)))
		end select

	end subroutine
