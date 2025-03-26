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

	subroutine process_file(filename, ttl, actl, string)

		! outputs: 
		! string: Array which includes all input descriptions
		! actl: Actual lines without any comments or tabs

		! inputs:
		! filename: File name 
		! ttl :: Total lines

		implicit none

		integer :: line, i, j
		integer, intent(in) :: ttl
		integer, intent(out) :: actl
		character(len=*), intent(in) :: filename
		character(len=200), intent(out) :: string(ttl)
		character(len=200) :: dummy_string(ttl)
		character(len=200) :: tag
		character, parameter :: htab = achar(9)


		open(1, file=filename, status='old', action='read')

		string = ''; dummy_string = ''
		do line = 1, ttl

			read(1, '(a)') string(line)
			tag = string(line)
			scribbing : do i = 1, len(tag)
				if(tag(i:i).eq.htab) then
					tag(i:i)= ''
				end if
				if(tag(i:i).eq.'!'.or. tag(i:i).eq.'#') then
					string(line) = tag(:i-1)
					exit scribbing
				end if
				string(line) = tag(:i)
			end do scribbing

		end do
		close(1)

		actl = 0
		do i = 1, ttl
			if(string(i).ne.'') then
				actl = actl + 1
				dummy_string(actl) = string(i)
			end if
		end do

		!dummy_string :: processed tags
		string = dummy_string

	end subroutine process_file