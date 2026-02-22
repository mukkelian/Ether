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

        subroutine read_input

	use init

	implicit none

	integer :: i, j, k, total_lines, actual_lines, txti
	logical :: file_present
	character(len=2) :: ab
	character(len=200) :: text, key, value

	call defaults	! Load default values

	inquire(file='in.ether', exist=file_present)
	if(.not.file_present) then
	if (rank == 0) then
		write(6, *) "==> Required file 'in.ether' is not present"
		write(6, *) "	 STOPPING now"
		write(6, *) ""
	end if
		stop
	end if

        open(unit=0, file='in.ether', status='old', action='read')

        total_lines = 0
        get_total_lines: do
        	read(0, "(a)", end=1) text
        	total_lines = total_lines + 1
        end do get_total_lines
1	close(0)

	allocate(input_data(total_lines))
	call process_file('in.ether', total_lines, actual_lines, input_data)

	do i = 1, actual_lines
		text = input_data(i)
		call lu(text, text, "L")
		j = index(text, '=')
		if(j.le.0) call at_not_equality(text)
		key = text(:j-1)
		value = text(j+1:)
!print*, "Key: ", trim(adjustl(key)), ", value: ", trim(adjustl(value))
		select case(trim(adjustl(key)))
		case("model")
			read(value, *) model
			Ising = .FALSE.; XYZ = .FALSE.
			if(model.eq.'ising') Ising = .TRUE.
			if(model.eq.'xyz') XYZ = .TRUE.
		case("mcs")
			read(value, *) tmcs, tmcs_eq
		case("temp")
			read(value, *) ht, lt, tint
		case("spin")
			read(value, *) (s(k), k = 1, nspecies)
		case("species")
			call count_species(value, n_speci_incl)
			allocate(species_to_include(n_speci_incl), stgg(n_speci_incl))
			read(value, *) (species_to_include(k), k = 1, n_speci_incl)
			do k = 1, n_speci_incl
				ab = species_to_include(k)
				call lu(ab(1:1), ab(1:1), 'U' )
				call lu(ab(2:2), ab(2:2), 'L' )
				species_to_include(k) = ab
			end do
		case("sc")
			read(value, *) sc(1:3)
			if(mod(sc(1), 2).ne.0.or.&
				mod(sc(2), 2).ne.0.or.&
				mod(sc(3), 2).ne.0) then
				write(6, *) ''
				write(6, *) '	ERROR'
				write(6, *) '	~~~~~'
				write(6, *) ''
				write(6, *) '	Dimension of supercell along x, y, and z &
				must be in the multiples of 2'
				write(6, *) '	STOPPING now'
				write(6, *) ''
				stop
			end if
		case("stg")
			read(value, *) staggered
		case("bc")
			read(value, *) bc(1:3)
		case("repeat")
			read(value, *) repeat
		case("angle")
			read(value, *) angle, dphi
			dphi = convert_to_rad*dphi
		case("coa")
			read(value, *) to_cal
		case("zeeman")
			read(value, *) Zeeman, h(1:3)
		case("g_factor")
			read(value, *) g_factor
		case("sia")
			read(value, *) sia
		case("para")
			read(value, *) para, J_para
		if((J_para.eq.0).and.para.and.rank == 0) then
			write(6, *) ""
			write(6, *) "ERROR:: parameter value is set ON and value &
       				cannot be zero. Kindly have a look into the input file"
			write(6, *) ""
			STOP
		end if
		case("ovrr")
			read(value, *) ovrr, ovrr_para, ovrr_MCS
			ovrr_para = abs(ovrr_para)
		case("seed")
			read(value, *) seed_value
			seed = seed_value

		case default
				print*, ''
				print*, ">	Found unknown/missing information in the given line"
				print*, "	      ~~~~~~~~~~~~~~~"
				print*, "	>> '",trim(adjustl(text)), "'"
				print*, ''
				print*, '	STOPPING now'
				print*, ''
				stop
			
		end select

	end do

	contains

	subroutine at_not_equality(remark)

		character(len=*), intent(in) :: remark
		print*, ''
		print*, ">	Missing '=' in the given line'"
		print*, "line: ", trim(adjustl(remark))
		print*, "	STOPPING now"
		print*, ""
		stop

	end subroutine at_not_equality

        end subroutine read_input
