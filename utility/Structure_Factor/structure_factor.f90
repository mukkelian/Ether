! Code Ether, based on the Monte Carlo technique, can be used to
! study the static and dynamics of spin models applied to any lattice geometry.
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

	implicit none
        integer :: nscan, sc(3), lattice_per_unit_cell, count, &
		i, j, k, l, ii, jj, kk, ll, nscan_qvec, iscan_qvec, total_ions, &
		p, p_, n_species, nbd_cell_x, nbd_cell_y, nbd_cell_z, fromx, &
        	fromy, fromz, tox, toy, toz

        real :: dx, dy, dz, lp1(3), lp2(3), vector(3), deg, &
		qvec, qvec_f, qvec_i, qvec_int, temp, dis(3), a(3, 3), &
        	to_angle, pi = 4.0*atan(1.0), light(3), mod_light, mod_vector, &
		theta, q_vec, r(3), k_vec(3), s_i(3,1), s_j(3,1), volume, &
		b(3,3), scaling(3), h_t, l_t, t_int

        real, allocatable :: ion(:, :, :, :, :)
	integer, allocatable :: site(:, :), ionn(:)
	complex :: img = (0,1), c_q
	character (len=80) :: lbl
	logical :: file_present = .FALSE.
        DOUBLE PRECISION :: start, finish
        call cpu_time(start)

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
        open(unit=1, file='piller.ether', status='unknown')
        open(unit=2,file='outputSF.dat', status='unknown',&
        action='write')
	open(unit=3, file='plot_sf.sh', status='unknown', action='write')
        open(unit=4, file='gss.dat', status='old', action='read')

        read(4, *) nscan, sc
        read(4, *) nbd_cell_x, nbd_cell_y, nbd_cell_z
        read(4, *) fromx, fromy, fromz, tox, toy, toz		
	read(4, *) lattice_per_unit_cell, n_species, h_t, l_t, t_int
	read(4, *) lbl

        allocate(ionn(n_species))
        read(4, *) (ionn(i), i = 1, n_species)
        structure : do i = 1,3
        	read(4,*) (a(i,j), j = 1,3) !lattice vectors
        end do structure
        allocate(ion(0:6, sc(1)+2*nbd_cell_x, sc(2)+2*nbd_cell_y, sc(3)+2*nbd_cell_z, lattice_per_unit_cell) )

        read(0, *) dx, dy, dz   ! directions
        read(0, *) deg
        read(0, *) qvec_f, qvec_i, qvec_int ! q vector
        read(0, *) scaling
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

        to_angle = 180./pi
        light = (/ dx, dy, dz/)

        ! READING LATTICES
        read(4,*)
	do l = 1, lattice_per_unit_cell
		do k = 1, sc(3) + 2*nbd_cell_z
			do j = 1, sc(2) + 2*nbd_cell_y
				do i = 1, sc(1) + 2*nbd_cell_x

                                        read(4,*) ion(4:6, i, j, k, l), ion(0, i, j, k, l)

                                end do
                        end do
                end do
        end do

	print*, 'DONE, now evaluatng structure factors (SF) as per given inputs'

        print*, ''
        print*, '::::::::::::::: STRUCTURE FACTOR (SF) :::::::::::::::'
        print*, ''
        print"(' SF along ',3f7.3 )",dx, dy,dz
        print"(' with scaling factors (x, y, z) :',3f7.3 )",scaling

        !SEARCHING LATTICE POINTS ALONG GIVEN DIRECTION
	count = 0; p_ = 0
	do l = 1, lattice_per_unit_cell
		do k = fromz, toz
			do j = fromy, toy
			        do i = fromx, tox

        lp1(1:3) = ion(4:6, i, j, k, l)
	p = 0
	do ll = 1, lattice_per_unit_cell
		do kk = fromz, toz
			do jj = fromy, toy
			        do ii = fromx, tox

	        lp2(1:3) = ion(4:6, ii, jj, kk, ll)

	        vector = lp1 - lp2
	
	        mod_light = sqrt(dot_product(light, light))
	        mod_vector = sqrt(dot_product(vector, vector))
	
	        theta = acos(sum(light*vector)/(mod_light*mod_vector))*to_angle

	        if((theta .le. deg).or.((180-theta) .le. deg))then
			p = 1
                          write(1, '(8I4, 4X, 2I10)') i, j, k, l, ii, jj, kk, ll,&
                          int(ion(0, i, j, k, l)), int(ion(0, ii, jj, kk, ll))
			count = count + 1
		end if

				end do
			end do
		end do
	end do

	if(p.eq.1) p_ = p_ + 1	!for counting same ion within piller

				end do
			end do
		end do
	end do

	close (1)

	open(unit=1, file='piller.ether', status='old', action='read')
	
	if(count==0) then
		print*, ''
		print*, 'ATTENTION!! no piller found'
		print*, 'STOPPING NOW'
		print*,''
		stop
	end if

	allocate(site(count,8))

	do i = 1, count
		read(1, *) site(i, :)
	end do
	close (1)
	call system('rm -rf piller*')

        !STRUCTURE FACTOR CORRELATION
        !################################################################

        nscan_qvec = abs(nint((qvec_f - qvec_i)/(qvec_int))) + 1

        reading_files : do kk = 1, nscan

        read(4, *) temp

       	do l = 1, lattice_per_unit_cell
       		do k = fromz, toz
       			do j = fromy, toy
				do i = fromx, tox
				
					read(4,*) ion(1:3, i, j, k, l)

				end do
			end do
		end do
	end do

	q_vector: do iscan_qvec = 1, nscan_qvec
		c_q = 0.0
	        q_vec = qvec_i + (qvec_int*(iscan_qvec-1)) ! for lower to higher

	do i = 1, count

		r = &
		ion(4:6, site(i, 1), site(i, 2), site(i, 3), site(i, 4)) - &
		ion(4:6, site(i, 5), site(i, 6), site(i, 7), site(i, 8))

		c_q = c_q + dot_product(&
		ion(1:3, site(i, 1), site(i, 2), site(i, 3), site(i, 4)), &
		ion(1:3, site(i, 5), site(i, 6), site(i, 7), site(i, 8)) &
		)* &
		exp(img*q_vec* dot_product(k_vec, r))

	end do

        !write(2,*)temp, q_vec, (abs(c_q) + p_)/(1.0*count+ p_) !p_ is the total number of sites
        write(2,*)temp, q_vec, (abs(c_q) + p_)/(1.0*p_) !p_ is the total number of sites

        end do q_vector

        write(2,*) ''

	end do reading_files

	close (2)

        !################################################################

	call cpu_time(finish)
	print*, ''
        print*, 'Creating polt file'
	write(3, *) '#set parameters accordingly'
	write(3, *) 'set nokey'
	write(3, *) 'set terminal png size 1100, 900 font "Times-New-Roman,22"'
	write(3, *) "set output 'figSF.png'"
	write(3, *) 'set tics font "Times-New-Roman,20"'
	write(3, *) 'set format y "%g"'
        write(3, 100)  l_t-t_int, abs(h_t-l_t+t_int)/5., h_t
	write(3, *) 'set mxtics 10'
        write(3, 101) qvec_f
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
        print*,"May the force be with you"
        print*,"~ Mukesh Kumar Sharma"
        print*,'e-mail: msharma1@ph.iitr.ac.in'
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
