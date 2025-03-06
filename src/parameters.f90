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

	subroutine parameters

	use init

	implicit none

	if (para) then
		lbl = ' KbT/J'
		j_exc = j_exc/para_value
		if(single_ion_anisotropy.and..not.Ising) sia_vec = sia_vec/para_value
		kb = real(1.0, dp)
		mb = real(1.0, dp)
		g_factor = real(1.0, dp)

		if(Zeeman) h = h/para_value
	else
		lbl = " K"
		! Converting 'j_exc' from meV into eV
		j_exc = j_exc/real(1000.0, dp)
		if(single_ion_anisotropy.and..not.Ising) sia_vec = sia_vec/real(1000.0, dp)
		! Bhor magneton
		mb = 5.7883818060d-5
		if(rank == 0) then
                write(6, *) ''
                write(6, *) "REMEMBER! 	For 'J in meV' case."
                write(6, *) '		ALL variables should be provided in the terms of milli orders (meV), '
                write(6, *) "		in 'input' file. eg., for 1meV put only 1 not 0.001"
                write(6, *) ''
                write(6, *) '		*******************************'
                write(6, *) "		KINDLY CHECK THE 'input' FILE"
                write(6, *) '		(IGNOR, if it is already DONE!)'
                write(6, *) '		*******************************'
                write(6, *) ''
                end if
	end if

	end subroutine parameters
