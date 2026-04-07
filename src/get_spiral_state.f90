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

	subroutine get_spiral_state(beta_value, observable_case)

	use init
	use omp_lib
		
	implicit none
	
	integer :: i, j, found_ion, dim1(2), dim2(2), &
		ions_for_spiral(spiral_capacity)

	real(dp), intent(in) :: beta_value
	real(dp) :: U_spiral_state, spiral_state_per_site, &
		vector1(3), vector2(3), uni_vec(3), ss(3), dis, r21(3), &
		spiral_state, spiral_state_chi, cp(3), ss_angle
		
	character(len=*), intent(in) :: observable_case
	character (len=200) :: remark
	
	complex(dp) :: ss_local(3)
	
	dim1 = (/1, spiral_capacity/)
	dim2 = (/1, 1/)
		
	select case(observable_case)
		
	case('ions_for_spiral')
	
	if(allocated(ifs))deallocate(ifs)
	allocate(ifs(spiral_capacity, total_ions))
	
	!$OMP PARALLEL DEFAULT(shared) PRIVATE(i,found_ion,vector1,j,vector2,&
	!$OMP dis,r21,uni_vec,ions_for_spiral)
	ions_for_spiral = 0
	!$OMP DO SCHEDULE(DYNAMIC)
	site_i: do i = 1, total_ions
		found_ion = 0
		vector1 = ion(6:8, i)
		site_j: do j = 1, total_ions
		if(j.ne.i) then
			vector2 = ion(6:8, j)

			r21 = vector2 - vector1
			dis = sqrt(dot_product(r21, r21))
			if(dis .le. ss_dis) then
				uni_vec = r21 / dis
				ss_angle = acos(dot_product(uni_vec, ss_direc))*to_angle
				if(ss_angle .le. ss_latency) then
				found_ion = found_ion + 1
				ions_for_spiral(found_ion) = j
				end if
			end if
		 end if
		 end do site_j

		ifs(:, i) = ions_for_spiral
	end do site_i
	!$OMP END DO
	!$OMP END PARALLEL
			
	case('spiral_state')

		!$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j,vector1,vector2,&
		!$OMP ions_for_spiral, cp) &
		!$OMP REDUCTION(+:ss)
		! Spiral State (SS)
		ss = 0.0_dp
		!$OMP DO SCHEDULE(DYNAMIC)
		calculate_ss: do i = 1, total_ions

			ions_for_spiral = ifs(:, i)
			vector1 = ion(1:3, i)

			get_ss_state: do j = 1, spiral_capacity
				if (ions_for_spiral(j) /= 0) then
					vector2 = ion(1:3, ions_for_spiral(j))
					call cross_product(vector1, vector2, cp)
					ss = ss + cp
				else
					exit get_ss_state
				end if
			end do get_ss_state

		end do calculate_ss
		!$OMP END DO

		!$OMP END PARALLEL

		! Get spiral state w.r.t. projection vector ss_proj
		spiral_state = sqrt(abs(dot_product(ss, ss_proj)))
			
		spiral_state_per_site = spiral_state/total_ions
		spiral_state_avg = spiral_state_avg + spiral_state_per_site
		spiral_state2_avg = spiral_state2_avg + spiral_state_per_site**2
		spiral_state4_avg = spiral_state4_avg + spiral_state_per_site**4

	case('avg_spiral_state')

		spiral_state_avg = spiral_state_avg/total_calculations
		spiral_state2_avg = spiral_state2_avg/total_calculations
		spiral_state4_avg = spiral_state4_avg/total_calculations

		s_spiral_state_avg = spiral_state_avg + s_spiral_state_avg
		e_spiral_state_avg(repeati) = spiral_state_avg

		U_spiral_state = 1.0_dp - (1.0_dp/3.0_dp)*(spiral_state4_avg/(spiral_state2_avg**2))
		s_U_spiral_state = U_spiral_state + s_U_spiral_state
		e_U_spiral_state(repeati) = U_spiral_state

		spiral_state_chi = (beta_value)*(spiral_state2_avg - spiral_state_avg**2)*total_ions
		s_spiral_state_chi = spiral_state_chi + s_spiral_state_chi
		e_spiral_state_chi(repeati) = spiral_state_chi

	case default

		remark = "Found unknown case tag '"//observable_case//&
		"' in get_spiral_state"
		call terminate (remark)

	end select

	end subroutine get_spiral_state
