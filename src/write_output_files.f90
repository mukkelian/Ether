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

	subroutine write_output_files(tempi)

		use init

		implicit none
		
		integer :: i
		integer, intent(in) :: tempi

		if(tempi.lt.1) then

			open(unit=20001, file='magnetization.dat', status='replace',action='write')
			open(unit=20002, file='energy.dat', status='replace', action='write')
			open(unit=20003, file='moment.dat', status='replace', action='write')

			m_head(1) ='# 1 Temp.'
			m_head(2) ='2 |M|'; m_head(3)='3 χ'
			m_head(4) ='4 Δ|M|'; m_head(5) ='5 Δχ'
			m_head(6) ='6 U(|M|)'; m_head(7) ='7 ΔU(|M|)'
		
		        e_head(1) ='# 1 Temp.'
			e_head(2) ='2 E'

	        	if(para)then
	        		e_head(3)='3 Cv' ; e_head(5) ='5 ΔCv'
			else
	                	e_head(3)='3 Cv/kB' ; e_head(5) ='5 Δ(Cv/kB)'
	        	end if

			e_head(4) ='4 ΔE'
			e_head(6) ='6 U(E)'; e_head(7) ='7 ΔU(E)'
			e_head(8) ='8 Accpt. (%)'

			write(20001,103) (m_head(i), i = 1, 7) 
			write(20002,104) (e_head(i), i = 1, 8)
103   			format(A13,A13,A13,1x,A13,1x,A13,2x,A13,A13)
104   			format(A13,A13,A13,A13,1x,A13,1x,A13,A13,A13)

			write(20003, *) '# spN = Nth species, M = magnitude of magnetic moment vector (Mx, My, Mz)'
			write(20003, *) '# Temp, sp1(Mx,My,Mz), sp1(Mx,My,Mz),..spN(Mx,My,Mz), sp1|M|, sp2|M|,..spN|M|'

			return
		end if

		!MAGNETIC MOMENT
        	write(20003, 101) temp_T(tempi), (mm_vector_avg_T(i, 1:3, tempi)&
        		/(product(sc)*ionn(i)), i = 1, nspecies), &
        	(sqrt(dot_product(mm_vector_avg_T(i, 1:3, tempi), mm_vector_avg_T(i, 1:3, tempi)))&
        		/(product(sc)*ionn(i)), i = 1, nspecies)
       
		!MAGNETIC
		write(20001,102) temp_T(tempi), s_mag_avg_T(tempi), &
		s_chi_T(tempi), err_mag_avg_T(tempi), &
		err_chi_T(tempi), s_U_mag_T(tempi), err_U_mag_T(tempi)

        	!ENERGY
        	write(20002,102) temp_T(tempi), s_eng_avg_T(tempi), &
		s_cv_T(tempi), err_eng_avg_T(tempi), &
		err_cv_T(tempi), s_U_eng_T(tempi), err_U_eng_T(tempi), acceptance_ratio(tempi)

102   		format(f11.5, 1x, es12.5, 1x,es12.5, 1x,es12.5, 1x,es12.5, 1x,es12.5, 1x,es12.5, 1x, f7.2)       
101		format(f11.5, 2X,500f9.3)

	end subroutine write_output_files
