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


	subroutine frac_to_cart(lattice_matrix, frac_coords, cart_coords)
	
	! Convert fractional coordinates to Cartesian coordinates using matrix multiplication
	! 
	! Arguments:
	!   lattice_matrix - 3x3 matrix with lattice vectors as columns [a b c]
	!                   Format: [[a_x, b_x, c_x], [a_y, b_y, c_y], [a_z, b_z, c_z]]
	!   frac_coords    - 3-element array [x, y, z] in fractional coordinates
	!   cart_coords    - 3-element array [X, Y, Z] in Cartesian coordinates (output)
	!
	! Formula: r_cart = L * r_frac
	
	use init

	real(dp), intent(in)  :: lattice_matrix(3, 3)  ! 3x3 lattice vectors matrix
	real(dp), intent(in)  :: frac_coords(3)       ! Fractional coordinates [x, y, z]
	real(dp), intent(out) :: cart_coords(3)       ! Cartesian coordinates [X, Y, Z]
	real(dp) :: r_frac(3)
	real(dp) :: r_cart(3)
	integer :: i

	r_frac = frac_coords

	! Performing matrix multiplication: r_cart = L * r_frac    
	do i = 1, 3
		r_cart(i) = lattice_matrix(1,i) * r_frac(1) + &
		lattice_matrix(2,i) * r_frac(2) + &
		lattice_matrix(3,i) * r_frac(3)
	end do

	cart_coords = r_cart

	end subroutine frac_to_cart
