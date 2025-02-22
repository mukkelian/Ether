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

      program vector
        implicit none
        integer :: i, j, num_ion, total_lines
        character (len=2), allocatable :: ion(:)
        real, allocatable :: a(:, :)
        character (len=2) :: ion_lbl, out, dummy
        character (len=2), allocatable :: ion_ID(:),  continue_tag(:)
        character(len=6) :: file_name
        logical :: file_present

        call system('cp *spK* sp_vector.xsf')

        inquire(file='sp_vector.xsf', exist=file_present)
	if(.not.file_present) then
		write(6, *) "==> spin file '*spK*' is absent"
		write(6, *) "	 STOPPING now"
		write(6, *) ""
		stop
        end if

        open (1, file='sp_vector.xsf', status='old', action='read')
        open (2, file='input', status='old', action='read')
        open (4, file='plot_vector.sh', status='unknown', action='write')

        read(2, *) num_ion ! number of ions
        allocate(ion(num_ion))
        read(2, *) (ion(i), i = 1, num_ion)
        close(2)

        do i = 1, num_ion
                dummy = ion(i)
                call lu(dummy(1:1), out, 'U' )  !U: upper case
		ion_lbl(1:1) = out
                call lu(dummy(2:2), out, 'L' )  !L: lower case
		ion_lbl(2:2) = out
                ion(i) = ion_lbl
        end do

        do i = 1, 10
                read(1, *)
        end do

        read(1, *) total_lines, i

        allocate(ion_ID(total_lines), a(total_lines, 6), continue_tag(num_ion))
        
        file_read : do i = 1, total_lines
                read(1, *) ion_ID(i), (a(i, j), j = 1, 6)
        end do file_read
        close(1)


        do i = 1, num_ion
                file_name = trim(adjustl(ion(i)))//'_vec'
                open(20+i, file=file_name, status='unknown', action='write')
                write(20+i, *) "# vector's data for given ions:"
                write(20+i, '(" #",5A3)') ion(i)
                do j = 1, total_lines
                if(ion(i).eq.ion_ID(j)) then
                        write(20+i, '(f7.4,2x,f7.4,2x,f7.4,2x,f7.4,2x,f7.4,2x,f7.4)') &
                                0.0, 0.0, 0.0, a(j, 4:6)
                end if
                end do
                close(20+i)
        end do

        print*, ''
        print*, 'DONE'
        print*, "Kindly check the gnuplot script, 'plot_vector.sh', &
                file"
        print*, "to plot the vectors"

        write(4, *) 'set parametric'
        write(4, *) 'set isosample 25,12'
        write(4, *) 'set ticslevel 0'
        write(4, *) 'R = 1.   # radius of sphere'
	write(4, *) 'set angle degree'
        write(4, *) 'set urange [0:360]'
        write(4, *) 'set vrange [-90:90]'
        write(4, *) 'set xlabel "X"'
        write(4, *) 'set ylabel "Y"'
        write(4, *) 'set zlabel "Z"'
        write(4, *) 'splot R*cos(u)*cos(v),R*cos(u)*sin(v),&
                     R*sin(u) w l t "" lc rgb "grey", \'
        continue_tag = ',\'
        continue_tag(num_ion) = ''
        do i = 1, num_ion
        write(4, "(' ',A8,' u 1:2:3:4:5:6 with vectors head filled lw 0.5 lc ',i3,1x,'t ', A10,1x,2A)") &
                '"'//trim(adjustl(ion(i)))//'_vec"', i, &
                '"'//trim(adjustl(ion(i)))//'_{vec}"', continue_tag(i)
        end do
        close(4)
        write(*, *) 'Or, '
        write(*, *) ''
        write(*, *) 'Use the given script'
        write(*, *) '~~~~~~~~~~~~~~~~~~~~'
        write(6, *) ''
        call system ('more plot_vector.sh')
        write(*, *) ''
        write(*, *) "May the force be with you"
        write(*, *) "~ Mukesh Kumar Sharma"
        write(*, *) "e-mail@ msharma1@ph.iitr.ac.in"
        write(*, *) ''

        contains

        !SMALL CAPS <----> LARGE CAPS #########################################
	!######################################################################

	subroutine lu(txtR, txtL, case_LU)
		character (len=*), intent(in) :: txtR
		character (len(txtR)), intent(out) :: txtL
        	character(len=53):: s
		character, intent(in) :: case_LU
		integer :: i, j
        	s = ' abcdefghijklmnopqrstuvwxyz ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		txtL = txtR
        	do i = 1, len(txtR)
        	        do j = 1, len(s)
        	                if( txtR(i:i).eq.s(j:j) )then
        	                        if( (j.le.27) .and. (case_LU.eq.'U') ) txtL(i:i) = s(j+27:j+27) ! into upper case
        	                        if( (j.gt.27) .and. (case_LU.eq.'L') ) txtL(i:i) = s(j-27:j-27) ! into lower case
        	                end if
        	        end do
        	end do
	end subroutine lu

	!######################################################################

      end program vector
