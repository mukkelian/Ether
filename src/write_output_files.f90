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

	subroutine write_output_files(T)

	use init

	implicit none
		
	integer :: i
	integer, intent(in) :: T

	if(T.lt.1) then

		open(unit=20001, file='magnetization.dat', status='replace',action='write')
		open(unit=20002, file='energy.dat', status='replace', action='write')
		open(unit=20003, file='moment_vectors.dat', status='replace', action='write')
		open(unit=20004, file='spiral_state_parameters.dat', status='replace', action='write')
		
		m_head(1) ='# 1 Temp.'
		m_head(2) ='2 |M|'; m_head(3)='3 χ'
		m_head(4) ='4 Δ|M|'; m_head(5) ='5 Δχ'
		m_head(6) ='6 U(|M|)'; m_head(7) ='7 ΔU(|M|)'
		
	        e_head(1) ='# 1 Temp.'
		e_head(2) ='2 E'
		e_head(4) ='4 ΔE'
		e_head(6) ='6 U(E)'; e_head(7) ='7 ΔU(E)'
		e_head(8) ='8 Accpt. (%)'

		if(ssp.and..not.ISING) then
		ss_head(1) ='# 1 Temp.'
		ss_head(2) ='2 |ψ|'; ss_head(3)='3 ψ'
		ss_head(4) ='4 Δ|ψ|'; ss_head(5) ='5 Δψ'
		ss_head(6) ='6 U(|ψ|)'; ss_head(7) ='7 ΔU(|ψ|)'
		end if

        	if(para)then
        		e_head(3)='3 Cv' ; e_head(5) ='5 ΔCv'
		else
                	e_head(3)='3 Cv/kB' ; e_head(5) ='5 Δ(Cv/kB)'
        	end if

		write(20001,103) (m_head(i), i = 1, 7) 
		write(20002,104) (e_head(i), i = 1, 8)
		if(ssp.and..not.ISING) write(20004,103) (ss_head(i), i = 1, 7)

		write(20003, *) '# Showing moment vectors (Mx, My, Mz) for species(sp):'
		write(20003, 105) (species(included_species_ID(i)), i = 1, total_species_to_include)
		write(20003, *) '# Temp, sp1(Mx,My,Mz), sp1(Mx,My,Mz),..spN(Mx,My,Mz), sp1|M|, sp2|M|,..spN|M|'

103   		format(A13,A13,A13,1x,A13,1x,A13,2x,A13,A13)
104   		format(A13,A13,A13,A13,1x,A13,1x,A13,A13,A13)
105		format(' #	', 5A3)
		return
	end if

	!MAGNETIC MOMENT VECTORS
       	write(20003, 101) temp_T(T), (mm_vector_avg_T(included_species_ID(i), 1:3, T)&
       		/tions(included_species_ID(i)), i = 1, total_species_to_include), &
       	(sqrt(dot_product(mm_vector_avg_T(included_species_ID(i), 1:3, T), &
                        mm_vector_avg_T(included_species_ID(i), 1:3, T))) &
       		/tions(included_species_ID(i)), i = 1, total_species_to_include)
       
	!MAGNETIC
	write(20001,102) temp_T(T), s_mag_avg_T(T), &
			s_chi_T(T), err_mag_avg_T(T), &
			err_chi_T(T), s_U_mag_T(T), err_U_mag_T(T)

	!ENERGY
       	write(20002,102) temp_T(T), s_eng_avg_T(T), &
			s_cv_T(T), err_eng_avg_T(T), &
			err_cv_T(T), s_U_eng_T(T), err_U_eng_T(T), acceptance_ratio(T)

	! SPIRAL SPIN STATES (SSS)
	if(ssp.and..not.ISING) then
	write(20004,102) temp_T(T), s_spiral_state_avg_T(T), &
			s_spiral_state_chi_T(T), err_spiral_state_avg_T(T), &
			err_spiral_state_chi_T(T), s_U_spiral_state_T(T), err_U_spiral_state_T(T)
	end if

102	format(f11.5, 1x, es12.5, 1x,es12.5, 1x,es12.5, 1x,es12.5, 1x,es12.5, 1x,es12.5, 1x, f7.2)       
101	format(f11.5, 2X,500f9.3)

	end subroutine write_output_files
