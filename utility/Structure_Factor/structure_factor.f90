! Code Ether, based on the Monte Carlo technique, can be used to
! stuqy the static and qynamics of spin models applied to any lattice geometry.
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

        program structure_factor

	use omp_lib

	implicit none

        integer :: nscan, sc(3), lattice_per_unit_cell, &
		i, j, k, l, total_q, ithQ, &
		p, p_, n_species, fromx, fromy, fromz, tox, toy, toz, total_ions

        real :: qx, qy, qz, lp1(3), lp2(3), vector(3), deg, &
		qvec, q_stop, q_start, dq, temp, dis(3), a(3, 3), &
        	to_angle, pi = 4.0*atan(1.0), light(3), mod_light, mod_vector, &
		theta, q, r(3), k_vec(3), Sx, Sy, Sz, s_i(3,1), s_j(3,1), volume, &
		b(3,3), scaling(3), h_t, l_t, t_int

        real, allocatable :: ion(:, :)
	integer, allocatable :: ionn(:)
	complex :: img = (0,1), S_q_x, S_q_y, S_q_z, expo, Cq(3)
	character (len=80) :: lbl
	logical :: file_present = .FALSE.
        DOUBLE PRECISION :: start, finish

        start = omp_get_wtime()

	print*, ''
	inquire(file='in_SF.ether', exist=file_present)
	if(.not.file_present) then
		write(*, *) "==> 'in_SF.ether' is not present"
		write(*, *) "	  STOPPING now"
		write(*, *) ""
		stop
	end if

	if(file_present) print*, "'in_SF.ether' file is found, now reading..."
        open(unit=0, file='in_SF.ether', status='old', action='read')
        read(0, *) qx, qy, qz   ! q_vec directions
        read(0, *) deg
        read(0, *) q_stop, q_start, dq ! q vector
        read(0, *) scaling

        print*, ''
        print*, '::::::::::::::: STRUCTURE FACTOR (SF) :::::::::::::::'
        print*, ''
        print"(' SF along ',3f7.3 )",qx, qy,qz
        print"(' with scaling factors (x, y, z) :',3f7.3 )",scaling

	inquire(file='gss.dat', exist=file_present)        
        if(.not.file_present) then
        	print*, ''
        	print*, 'STOPPING now'
        	print*, "'gss.dat' file is not present"
        	print*, ''
        	stop
        end if

	print*, ''
        print*, "Calculating SF using given 'gss.dat' file"

	file_present = .FALSE.
        open(unit=4, file='gss.dat', status='old', action='read')

        read(4, *) nscan, sc
        read(4, *) fromx, fromy, fromz, tox, toy, toz		
	read(4, *) lattice_per_unit_cell, n_species, h_t, l_t, t_int
	read(4, *) lbl

        allocate(ionn(n_species))
        read(4, *) (ionn(i), i = 1, n_species)
        structure : do i = 1,3
        	read(4,*) (a(i,j), j = 1,3) !lattice vectors
        end do structure

        total_ions = product(sc)*lattice_per_unit_cell

        allocate(ion(0:6, total_ions))

        a(1, :) = a(1, :)*scaling(1)
        a(2, :) = a(2, :)*scaling(2)
        a(3, :) = a(3, :)*scaling(3)

        ! Evaluation of reciprocal vectors
        s_i = 0.0; s_j = 0.0
        s_i(:,1) = a(2,:); s_j(:,1) = a(3,:)

        volume = dot_product(a(1,:), cross_product(s_i(:,1),s_j(:,1)))

        s_i(:,1) = a(2,:); s_j(:,1) = a(3,:)
        b(1,:) =  cross_product(s_i(:,1),s_j(:,1))

        s_i(:,1) = a(3,:); s_j(:,1) = a(1,:)
        b(2,:) =  cross_product(s_i(:,1),s_j(:,1))

        s_i(:,1) = a(1,:); s_j(:,1) = a(2,:)
        b(3,:) =  cross_product(s_i(:,1),s_j(:,1))

        b = b*2.*pi/dble(volume)

        k_vec(:) = b(1,:) + b(2,:) + b(3,:)

        ! READING LATTICES
        read(4,*)
	do i = 1, total_ions
                read(4,*) ion(4:6, i), ion(0, i)
        end do

        !STRUCTURE FACTOR CORRELATION
        !################################################################

        total_q = abs(nint((q_stop - q_start)/(dq))) + 1

        open(unit=2,file='outputSF.dat', status='unknown', action='write')
        
        reading_files : do k = 1, nscan

        read(4, *) temp


	do i = 1, total_ions
		! Reading Spin vectors for i-th ion
		read(4,*) ion(1:3, i)
	end do

	q_vector: do ithQ = 1, total_q

	S_q_x = 0.0
	S_q_y = 0.0
	S_q_z = 0.0

	! for lower to higher
        q = q_start + (dq*(ithQ-1))
	
	!$omp parallel do default(shared) private(i, r, Sx, Sy, Sz, expo) &
	!$omp reduction(+:S_q_x, S_q_y, S_q_z)
	calculate_Sq: do i = 1, total_ions

		r(1:3) = ion(4:6, i)
		Sx = ion(1, i)
		Sy = ion(2, i)
		Sz = ion(3, i)

		expo = exp(img*q*dot_product(k_vec, r))

		S_q_x = S_q_x + &
		Sx*expo

		S_q_y = S_q_y + &
		Sy*expo

		S_q_z = S_q_z + &
		Sz*expo

	end do calculate_Sq
	!$omp end parallel do
	
	Cq = (/S_q_x, S_q_y, S_q_z/)

        write(2,*)temp, q, real(dot_product(Cq, Cq))/total_ions**2, &
        	real(abs(S_q_x))/total_ions, real(abs(S_q_y))/total_ions, &
        	real(abs(S_q_z))/total_ions

        end do q_vector

        write(2,*) ''

	end do reading_files

	close (2)

        !################################################################

	finish = omp_get_wtime()

	open(unit=3, file='plot_sf.sh', status='unknown', action='write')

	print*, 'Done'
	print*, ''
        print*, 'Creating polt file'
	write(3, *) '#set parameters accordingly'
	write(3, *) 'unset key'
	write(3, *) 'set terminal png size 1100, 900 font "Times-New-Roman,22"'
	write(3, *) "set output 'figSF.png'"
	write(3, *) 'set tics font "Times-New-Roman,20"'
	write(3, *) 'set format y "%g"'
        write(3, 100)  l_t-t_int, abs(h_t-l_t+t_int)/5., h_t
	write(3, *) 'set mxtics 10'
        write(3, 101) q_stop
	write(3, *) 'set mytics 5'
	write(3, *) "set title 'Structure Factor Vs. Temperature'"
	write(3, *) 'unset surface'
	write(3, *) 'unset key'
	write(3, *) 'set border lw 2'
	write(3, *) 'set palette defined (0 "blue", 0.4 "yellow", 1 "red")'
	write(3, *) 'set view 360.0, 359.99'
	write(3, *) 'set pm3d at b'
	if(trim(adjustl(lbl)).eq.'K')then
		write(3, 103) 'K'
	else
		write(3, 102) 'k_{B}T/J'
	end if
	write(3, *) "set ylabel 'Q vec.' tc 'red'"
	write(3, *) 'set xrange [:]'
	write(3, *) 'set yrange[:]'
	write(3, *) "#set zlabel 'C_{q}'"
	write(3, *) "splot 'outputSF.dat' u 1:2:3 w l lc 'grey' t ''"
	close(3)
	print*, ''
	print*, 'command to plot Structure Factor ==> "gnuplot plot_sf.sh"'
	print*, ''
        print*, 'Completed!'
	print*, ''
        print '(" Total calculation time = ",i4," Hr",i3," min.",i3,&
        " sec.")'&
        ,int((finish-start)/3600.0),int(mod((finish-start), &
        3600.0)/60.0),&
        int(mod(mod((finish-start), 3600.0),60.0))
        print*,''
        print*,''
100	format (" set xtics ", f11.5,",", f11.5,",", f11.5," scale &
        2, 0.75 textcolor 'blue'")
101	format (" set ytics 0,0.25,",f10.5," scale 1, 0 textcolor 'red'")
102	format (" set xlabel 'Temperature (",A8, ")' tc 'blue'")
103	format (" set xlabel 'Temperature (",A1, ")' tc 'blue'")
	contains

        FUNCTION cross_product(v1, v2)
 
 	       real, DIMENSION(3) :: cross_product
 	       real, DIMENSION(3), INTENT(IN) :: v1, v2

 	       cross_product(1) = (v1(2) * v2(3)) - v1(3) * v2(2)
 	       cross_product(2) = (v1(3) * v2(1)) - v1(1) * v2(3)
 	       cross_product(3) = (v1(1) * v2(2)) - v1(2) * v2(1)

        END FUNCTION cross_product

	end program structure_factor
