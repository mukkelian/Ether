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

	subroutine write_gss(tempi)
	
	use init
		
	implicit none
		
	integer :: i, j, k, l, m
	integer, intent(in) :: tempi
	gss_ID = 10009

	if(tempi.le.0) then

	open(unit=gss_ID, file='gss.dat', status='replace', action='write')

        if(para)then
                lbl = ' KbT/J'
        else
                lbl = " K"
        end if

	! STORING crystal information in GSS (Ground Spin States)
       	write(gss_ID, '(I5,1x, I5, 1x, I5, 1x, I5)') nscan, sc
       	write(gss_ID, '(I5,1x, I5, 1x, I5)') nbd_cell_x, nbd_cell_y, nbd_cell_z
       	write(gss_ID, '(I5,1x, I5, 1x, I5, 1x, I5, 1x, I5, 1x, I5)') fromx, &
       		fromy, fromz, tox, toy, toz
       	write(gss_ID, '(I5,1x, I5, f11.5, 1x, f11.5, 1x, f11.5,1x, A6)') &
       		lattice_per_unit_cell, nspecies, ht, lt, tint
       	write(gss_ID, '(A20)') lbl
       	write(gss_ID, '(10I8)') tions(1:)
	do i = 1,3
		write(gss_ID,'(3f10.5)') (abc(i,j), j= 1,3)
	enddo
	write(gss_ID,*) '# lattice points'
	do l = 1, lattice_per_unit_cell
		do k = 1, sc(3) + 2*nbd_cell_z
			do j = 1, sc(2) + 2*nbd_cell_y
				do i = 1, sc(1) + 2*nbd_cell_x

			write(gss_ID,'(f11.6, 1x, f11.6, 1x, f11.6, 1X, i8)') &
				ion(6:8, i, j, k, l), int(ion(0, i, j, k, l))

				end do
			end do
		end do
	end do
	return

	end if

        !GSS.dat
        m = 0
	write(gss_ID, '(f11.5)') temp_T(tempi)

	do l = 1, lattice_per_unit_cell
		do k = fromz, toz
			do j = fromy, toy
				do i = fromx, tox

			write(gss_ID, '(f11.6, 1x, f11.6, 1x, f11.6,1X, i8, 1x, i8)') &
			global_spn(m + 1 + lspn), global_spn(m + 2 + lspn), &
			global_spn(m + 3 + lspn), int(global_spn(m + 4 + lspn))
			m = m + 4

				end do
			end do
		end do
	end do

	if(tempi.eq.nscan) close(gss_ID)

	end subroutine write_gss
