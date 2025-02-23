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

	subroutine graph

		use init
		
		implicit none

	        if(para)then
	                lbl = ' k_{B}T/J'
	        else
	                lbl = "K"
	        end if

		open(10008, file='graph.sh', status='unknown')
		write(10008,*) '# set labels, lables, and xtics/ytics accordingly'
		write(10008,*) 'set nokey'
		write(10008,*) 'set terminal png size 1100, 900 font "Times-New-Roman,18"'
		write(10008,*) "set output 'results_Ether.png'"
		write(10008,*) 'set multiplot'                     
		write(10008,*) "set size squar 0.5,0.5"
		write(10008,*) 'set nokey'
		write(10008,*) 'set format y "%g"'
	        write(10008, 303)  lt - tint, abs(ht - lt + tint)/5., ht
		write(10008,*) 'set mxtics 10'
		write(10008,*) 'set origin 0.0,0.0'
		write(10008,*) 'set grid'
	        if(para)then
		        write(10008, 301) trim(lbl)
	        else
	                write(10008, 302) trim(lbl)
		end if
	        if(staggered) then
	        	write(10008,*) 'set ylabel "Staggered Magnetization (M)"'
	        else
	        	write(10008,*) 'set ylabel "Magnetization (M)"'
	        end if
	
		write(10008,*) "plot 'magnetization.dat' u 1:2:4 with yerrorbars &
	                        lt 6 notitle"
		write(10008,*) 'set origin 0.5, 0.0'
	        if(para)then
	                write(10008, 301) lbl
	        else
	                write(10008, 302) lbl
	        end if
	
	        if(staggered) then
	        	write(10008,*) 'set ylabel "Staggered Susceptibility (χ)"'
	        else
	        	write(10008,*) 'set ylabel "Susceptibility (χ)"'
	        end if
	        
		write(10008,*) "plot 'magnetization.dat' u 1:3:5 with yerrorbars &
	                        lt 6 notitle"
		write(10008,*) 'set origin 0.0,0.5'
	        if(para)then
	                write(10008, 301) lbl
	        else
	                write(10008, 302) lbl
	        end if
		write(10008,*) 'set ylabel "Energy (meV)"'
		write(10008,*) "plot 'energy.dat' u 1:($2*1000):4 with yerrorbars &
	                        lt 6 notitle"
		write(10008,*) 'set origin 0.5,0.5'
	        if(para)then
	                write(10008, 301) lbl
	        else
	                write(10008, 302) lbl
	        end if
	        if(para)then
	                write(10008, *) 'set ylabel "Specific Heat (C_{V})"'
	        else
	                write(10008, *) 'set ylabel "Specific Heat &
			(C_{V}k_{B}^{-1})"'
	        end if
		write(10008,*) "plot 'energy.dat' u 1:3:5 with yerrorbars &
	                        lt 6 notitle"
		write(10008,*) 'unset multiplot'
		close (10008)
	
301     	format(' set xlabel "Temperature (',A8,')"')
302     	format(' set xlabel "Temperature (',A1,')"')
303     	format (" set xtics ", f11.5,",", f11.5,",", f11.5)

	end subroutine graph
