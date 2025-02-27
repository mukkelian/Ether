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

	subroutine startup(values)

	use init

	implicit none

        integer,dimension(8), intent(in) :: values

	write(6, *) ""
	write(6, *) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"	        
	write(6, *) "~~~~~                      ETHER                        ~~~~~"
	write(6, *) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	write(6, *) ""
        write(6, '(" PROGRAM STARTED on date ",i2,"-",i2,"-",i4)') values(3),values(2), values(1)
	write(6, "(' at time ',i2,' hrs. ',i2,' min. ',i2,' sec. ')") values(5:7)
	write(6, *)""
	
	if(Ising) then
		write(6, *)"Based on Ising model (H = JSi.Sj)"
	elseif(Heisenberg) then
		write(6, *)"Based on classical Heisenberg model (H = JSi*Sj)"
	else
		write(6, *) "WARNING: unable to get right model ID"
		write(6, *) "setting model to default 'Heisenberg model (ID = H)'"
		model = 'h'
	end if
	
	write(6, *)""
	write(6, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	write(6, *) 'Working on Monte Carlo steps'
	write(6, *) 'per Spin (MCS)'
	write(6, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	write(6, "(1X,A21,I8)") 'Total MCS steps: ', tmcs
	write(6, "(1X,A21,I8)") 'Equilibration steps: ', tmcs_eq
	write(6, *) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	write(6, *) ''
        write(6, '(I4,1x,I4,1x,I4,"  :: Lx, Ly, Lz <==> Supercell size")') sc

        write(6, *)''
        write(6, '(2X,A2,3x,A2,3x,A2,2x"(boundary conditions along x,y, and z-axis, respectively)")') bc(1:3)
        write(6, *)"                    'o' => open; 'c' => closed"
        write(6, *)''
	write(6, *)'==> Spins are set in random configurations'

        if((angle.eqv..FALSE.).and.Heisenberg) then
        	write(6, *) '==> Spin vectors will be chosen as per George Marsaglia method'
        	write(6, *)''	
        	write(6, *) '     >  George Marsaglia, The Annals of Mathematical Statistics'
        	write(6, *) '        Vol. 43, No. 2, 645-646 (1972).'
        	write(6, *)''
        elseif(angle) then
        	write(6, '(" ==> +/-",f6.2," deg. angle has been chosen for the convergence")') dphi/convert_to_rad
        	write(6, *)''
        	write(6, *) '     >  Application of the Monte Carlo Method in Statistical Physics'
        	write(6, *) '        by Kurt Binder, Second Edition'
        	write(6, *)''
        end if

	if(Ising) then
		write(6, *) '==> Spin states will be chosen from two states i.e., +S and -S'
	end if
	
        if (ovrr.and..not.Ising)then
        	write(6, '(" ==> Overrelaxation method with ",I4," OVR steps has been considered along with &
        		Metropolis algorithm")') overrelaxed
        	write(6, *)''
        	write(6, *) '     >  Michael Creutz, Physical Review D Vol. 36, 2 (1987).'
        	write(6, *) '     >  J. L. Alonso and A. Tranacon, Physical Review B Vol. 53, 5 (1996).'
        	write(6, *) '     >  Nicholas Metropolis, Arianna W. Rosenbluth, Marshall N. Rosenblunth,'
        	write(6, *) '        Augusta H. Teller, and Edward Teller, The Journal of Chemical Physics'
        	write(6, *) '        Vol. 21, 1087 (1953).'
       		write(6, *)''
        else
        	write(6, *) '==> Metropolis algorithm will be used'
        	write(6, *)''
        	write(6, *) '     >  Nicholas Metropolis, Arianna W. Rosenbluth, Marshall N. Rosenblunth,'
        	write(6, *) '        Augusta H. Teller, and Edward Teller, The Journal of Chemical Physics'
        	write(6, *) '        Vol. 21, 1087 (1953).'
        	write(6, *)''
        end if

	if(EXalgo) then
		write(6, *) '==> Exchange alogorithm will be used along with Metropolis algorithm'
		write(6, '("     with exchange interval", i4)') exchange_interval
		write(6, *) ''
		write(6, *) '     >  Koji Hukushima and Koji Nemoto, Journal of the Physical Society'
		write(6, *) '	 of Japan 65, 1604-1608 (1996).'
		write(6, *) ''

	if(temp_ex.eqv..FALSE.)then
		write(6, '(" ==> Temperature points will be generated as per the method given in the appendix")')
		write(6, *) ''
		write(6, *) '     >  Koji Hukushima, Physical Review E 60, 4 (1999).'
		write(6, *)''
	else
		write(6, *) '==> Equi-spaced temperature points have been considered for Exchange algorithm'
	end if

	end if

	if(staggered)then
		write(6, *)'==> Staggered magnetization is SELECTED!'
	end if

	if(Zeeman)then
                if(para)then
                        lbl = ' as H/J'
                else
                        lbl = " in Tesla"
                end if
		write(6, '(" ==> Magnetic field ::", 3f10.5, " (Mx, My, Mz)",A9," is SELECTED!")') h, lbl
	end if

	if(single_ion_anisotropy)then
                if(para)then
                        lbl = " as SIA/J"
                else
                        lbl = " in meV"
                end if
		write(6, '(" ==> SIA (Single Ion Anisotropy) ", A9," is SELECTED!")') lbl
	end if

	if(para)then
                write(6, *) '==> Calculations in terms of parameters is ON!'
                write(6, '(" ==> Parameters w.r.t. ",f8.3, " meV has been taken")') para_value
	end if

        end subroutine startup
