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

	subroutine rw_file(mode, filename, file_unit, dim1, dim2, at, content)

	use init, only: dp

	implicit none

	integer(kind=8) :: total_record
	integer, intent(in) :: dim1(2), dim2(2), at, file_unit
	integer :: io_status, size

	real(dp), intent(inout):: content(dim1(1):dim1(2), dim2(1):dim2(2))

	character(len=*), intent(in) :: filename, mode

	logical :: found = .FALSE.

	inquire(iolength=total_record) content
	inquire(file=filename, exist=found, size=size)

	if(.not.found.and.mode.ne.'r') then
		open(file_unit, file=trim(adjustl(filename)), form='unformatted', &
		access='direct', recl=total_record, action='readwrite')
		close(file_unit)
	end if

	open(file_unit, file=trim(adjustl(filename)), form='unformatted', &
	access='direct', status='old', recl=total_record, iostat=io_status)

	if(mode.eq.'w') then
		write(file_unit, rec=at, iostat=io_status) content
		if (io_status /= 0) then
			print *, "	Error: ", io_status
			call terminate ('Error in writing data')
		end if
		close(file_unit)
	elseif(mode.eq.'r') then
	
		if(size.eq.0) then
			print *,"	'",trim(adjustl(filename)), "' file is empty"
		end if
		
		read(file_unit, rec=at, iostat=io_status) content
		if (io_status /= 0) then
			print *, "	Error: ", io_status
			call terminate ('Error in reading data')
		end if
		close(file_unit)
	else
		print*, ''
		print*, "Detected unknown input argument '", trim(adjustl(mode)), "' on rw_file subroutine"
		call terminate (mode)
	end if

	end subroutine rw_file
