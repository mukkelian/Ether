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

!=======================================================================
!  Reads the TB2J output file "exchange.out" and writes:
!
!     J_iso    : isotropic exchange J_iso      (scalar, meV)  - always
!     J_ani  : anisotropic exchange J_ani    (3x3,    meV)  - only if
!				"J_ani:" is present in exchange.out
!     DMI      : Dzyaloshinskii-Moriya vector  (3x1,    meV)  - only if
!				"DMI:" is present in exchange.out
!     J_prime  : J'						    (scalar, meV)  - only if
!				"Jprime:" is present in exchange.out
!
!  Arguments
!  ---------
!   bond_length (in) : real(8). Only ion pairs with distance
!				      <= bond_length (Angstrom) are written.
!   home		(in) : logical.
!				      .true.  -> only ion pairs with R = (0, 0, 0)
!								 (both ions in the home cell)
!				      .false. -> ion pairs from any cell
!=======================================================================

	subroutine extract_TB2J(bond_length, home)
	implicit none

	integer, parameter :: dp = kind(1.0d0)

	real(dp), intent(in) :: bond_length
	logical,  intent(in) :: home
	character(len=*), parameter :: infile = 'exchange.out'
	integer, parameter  :: unit_in = 41
	integer, parameter  :: maxlen  = 512
	real(dp), parameter :: tol_d   = 1.0d-4		! two distances are equal if |d1-d2| < tol_d
	real(dp), parameter :: tol_0   = 1.0d-6		! a vector component is zero if |v| < tol_0

	integer :: nbond	! number of bonds kept
	character(len=16), allocatable :: lab1(:), lab2(:)	! labels
	real(dp), allocatable :: vec(:,:)	! bond vector
	real(dp), allocatable :: dist(:)	! bond length
	real(dp), allocatable :: jiso(:)	! J_iso
	real(dp), allocatable :: jprime(:)	! J'
	real(dp), allocatable :: dmi(:,:)	! DMI vector
	real(dp), allocatable :: jani(:,:,:)	! J_ani
	integer,  allocatable :: bond_id(:)	! index of the bond length

	logical :: found_jprime, found_dmi, found_jani
	integer :: nshell
	real(dp), allocatable :: shell(:)
	character(len=maxlen) :: line
	logical  :: exists
	real(dp) :: cutoff
	integer  :: ios, nmax

	inquire(file=infile, exist=exists)
	if (.not. exists) then
		write(6,'(a)') '	> ERROR: the file "'//infile//'" is absent in the current directory.'
		stop 'extract_TB2J: exchange.out not found'
	end if

	if (bond_length <= 0.0_dp) then
		write(6,'(a)') '	> ERROR: bond_length must be a positive number (Angstrom).'
		stop 'extract_TB2J: wrong bond_length'
	end if
	cutoff = bond_length


	! count the separator lines "-----"
	open(unit=unit_in, file=infile, status='old', action='read')
	nmax = 0
	do
		read(unit_in,'(a)',iostat=ios) line
		if (ios /= 0) exit
		if (is_separator(line)) nmax = nmax + 1
	end do
	close(unit_in)
	nmax = max(nmax, 1)

	allocate(lab1(nmax), lab2(nmax), vec(3,nmax), dist(nmax), bond_id(nmax))
	allocate(jiso(nmax), jprime(nmax), dmi(3,nmax), jani(3,3,nmax))

	!read the bonds
	call read_bonds()

	if (nbond == 0) then
	write(6,'(a,f10.4,a)') '	> No bond found with distance <= ', cutoff, ' A.'
	end if


	!Distinct bond lengths (sorted) and bond_ID of every bond
	allocate(shell(max(nbond,1)))
	call find_shells()

	!Writing the output files
	if (home) then
		write(6,'(a)') '	> extract_TB2J: only ion pairs with R = (0, 0, 0)'
	else
		write(6,'(a)') '	> extract_TB2J: ion pairs with all R'
	end if
	write(6,'(a,i0,a,i0,a)') '	> extract_TB2J: ', nbond, ' bonds kept, ', nshell, &
		' distinct bond lengths.'

	call write_file('J_iso', 1)
	write(6,'(a)') '	> written : J_iso'

	if (found_jani) then
		call write_file('J_ani', 2)
		write(6,'(a)') '	> written : J_ani'
	else
		write(6,'(a)') '	> J_ani not present in exchange.out -> J_ani not created'
	end if

	if (found_dmi) then
		call write_file('DMI', 3)
		write(6,'(a)') '	> written : DMI'
	else
		write(6,'(a)') '	> DMI    not present in exchange.out -> DMI not created'
	end if

	if (found_jprime) then
		call write_file('J_prime', 4)
		write(6,'(a)') '	> written : J_prime'
	else
		write(6,'(a)') '	> Jprime not present in exchange.out -> J_prime not created'
	end if

	deallocate(lab1, lab2, vec, dist, bond_id, jiso, jprime, dmi, jani, shell)

	contains

	! Read exchange "Exchange:" block
	subroutine read_bonds()
	logical  :: in_exchange, after_sep, keep
	real(dp) :: d_tmp, j_tmp, v_tmp(3), m(3,3)
	integer  :: r_tmp(3)
	logical  :: in_range, in_home
	character(len=16) :: a1, a2
	character(len=maxlen) :: txt

	nbond = 0
	found_jprime = .false.
	found_dmi = .false.
	found_jani = .false.
	in_exchange = .false.
	after_sep = .false.
	keep = .false.

	open(unit=unit_in, file=infile, status='old', action='read')
	do
		read(unit_in,'(a)',iostat=ios) line
		if (ios /= 0) exit
		txt = adjustl(line)

		if (.not. in_exchange) then
			if (txt(1:9) == 'Exchange:') in_exchange = .true.
			cycle
		end if

		if (is_separator(line)) then
			keep = .false.
			after_sep = .true.
			cycle
		end if
		if (len_trim(txt) == 0) cycle

		!first line after a separator = bond line
		if (after_sep) then
			after_sep = .false.
			if (.not. is_bond_line(txt)) then
			if (index(txt,'(') == 0) exit
			cycle
		end if
		call parse_bond_line(txt, a1, a2, r_tmp, j_tmp, v_tmp, d_tmp)

		in_range = (d_tmp <= cutoff + tol_d)
		in_home  = all(r_tmp == 0)
		if (in_range .and. (in_home .or. .not. home)) then
			keep  = .true.
			nbond = nbond + 1
			lab1(nbond)     = a1
			lab2(nbond)     = a2
			vec(:,nbond)    = v_tmp
			dist(nbond)     = d_tmp
			jiso(nbond)     = j_tmp
			jprime(nbond)   = 0.0_dp
			dmi(:,nbond)    = 0.0_dp
			jani(:,:,nbond) = 0.0_dp
		end if
		cycle
		end if

		if (txt(1:6) == 'J_iso:') then
			if (keep) then
				read(txt(7:),*,iostat=ios) j_tmp
				if (ios == 0) jiso(nbond) = j_tmp
			end if

		else if (index(txt,'Jprime:') > 0) then
			found_jprime = .true.
			if (keep) jprime(nbond) = read_after(txt, 'Jprime:')

		else if (index(txt,'DMI:') > 0) then
		found_dmi = .true.
		call read_vector_in_brackets(txt, v_tmp)
		if (keep) dmi(:,nbond) = v_tmp

		else if (index(txt,'J_ani:') > 0) then
			found_jani = .true.
			call read_matrix(m)
			if (keep) jani(:,:,nbond) = m
		end if
	end do
	close(unit_in)
	end subroutine read_bonds

	subroutine read_matrix(m)
	real(dp), intent(out) :: m(3,3)
	character(len=maxlen) :: row
	integer :: k
	m = 0.0_dp
	do k = 1, 3
		read(unit_in,'(a)',iostat=ios) row
		if (ios /= 0) return
		call replace_chars(row, '[],', ' ')
		read(row,*,iostat=ios) m(k,1), m(k,2), m(k,3)
	end do
	end subroutine read_matrix

	subroutine read_vector_in_brackets(txt, v)
	character(len=*), intent(in) :: txt
	real(dp), intent(out) :: v(3)
	integer :: i1, i2
		v  = 0.0_dp
		i1 = index(txt,'(')
		i2 = index(txt,')')
		if (i1 == 0 .or. i2 <= i1) return
		read(txt(i1+1:i2-1),*,iostat=ios) v
	end subroutine read_vector_in_brackets

	function read_after(txt, key) result(x)
	character(len=*), intent(in) :: txt, key
	real(dp) :: x
	character(len=maxlen) :: rest
	integer :: i
	x = 0.0_dp
	i = index(txt, key)
	if (i == 0) return
	rest = txt(i+len(key):)
	i = index(rest, ',')
	if (i > 0) rest = rest(1:i-1)
	read(rest,*,iostat=ios) x
	if (ios /= 0) x = 0.0_dp
	end function read_after

	subroutine parse_bond_line(txt, a1, a2, r, j, v, d)
	character(len=*), intent(in)  :: txt
	character(len=16), intent(out) :: a1, a2
	integer, intent(out) :: r(3)
	real(dp), intent(out) :: j, v(3), d
	integer :: p1, q1, p2, q2

	p1 = index(txt, '(')				 ! first  "("  -> R
	q1 = index(txt, ')')
	p2 = p1 + index(txt(p1+1:), '(')     ! second "("  -> vector
	q2 = q1 + index(txt(q1+1:), ')')

	read(txt(1:p1-1),*)     a1, a2
	read(txt(p1+1:q1-1),*)  r
	read(txt(q1+1:p2-1),*)  j
	read(txt(p2+1:q2-1),*)  v
	read(txt(q2+1:),*)      d
	end subroutine parse_bond_line

	logical function is_bond_line(txt)
	character(len=*), intent(in) :: txt
	character :: c
	integer :: i1, i2
	is_bond_line = .false.
	c = txt(1:1)
	if (.not. ((c >= 'A' .and. c <= 'Z') .or. (c >= 'a' .and. c <= 'z'))) return
	if (txt(1:6) == 'J_iso:') return
 	i1 = index(txt, '(')
	if (i1 == 0) return
	i2 = index(txt(i1+1:), '(')
	if (i2 == 0) return
	is_bond_line = .true.
	end function is_bond_line

	logical function is_separator(txt)
	character(len=*), intent(in) :: txt
	is_separator = .false.
	if (len_trim(txt) < 10) return
	is_separator = (verify(trim(adjustl(txt)), '-') == 0)
	end function is_separator

	subroutine replace_chars(str, chars, by)
	character(len=*), intent(inout) :: str
	character(len=*), intent(in) :: chars
	character, intent(in)    :: by
	integer :: i
	do i = 1, len_trim(str)
		if (index(chars, str(i:i)) > 0) str(i:i) = by
	end do
	end subroutine replace_chars

	function element_name(label) result(name)
	character(len=*), intent(in) :: label
	character(len=16) :: name
	integer :: i, k
	name = ' '
	k = 0
	do i = 1, len_trim(label)
	if (label(i:i) < '0' .or. label(i:i) > '9') then
		k = k + 1
		name(k:k) = label(i:i)
	end if
	end do
	end function element_name

	integer function sign_of(x)
	real(dp), intent(in) :: x
	if (x > tol_0) then
		sign_of = 1
	else if (x < -tol_0) then
		sign_of = -1
	else
		sign_of = 0
	end if
	end function sign_of

	subroutine find_shells()
	integer :: ib, is, k
	logical :: found
	real(dp) :: x

	! distinct values
	nshell = 0
	do ib = 1, nbond
		found = .false.
		do is = 1, nshell
			if (abs(dist(ib) - shell(is)) < tol_d) then
				found = .true.
				exit
			end if
		end do
		if (.not. found) then
			nshell = nshell + 1
			shell(nshell) = dist(ib)
		end if
	end do

	! sort ascending
	do is = 2, nshell
		x = shell(is)
		k = is - 1
		do while (k >= 1)
			if (shell(k) <= x) exit
				shell(k+1) = shell(k)
				k = k - 1
		end do
		shell(k+1) = x
	end do

	! bond_ID of every bond
	do ib = 1, nbond
		bond_id(ib) = 0
		do is = 1, nshell
		if (abs(dist(ib) - shell(is)) < tol_d) then
			bond_id(ib) = is
			exit
		end if
		end do
	end do
	end subroutine find_shells

	subroutine write_file(fname, kind_of)
	character(len=*), intent(in) :: fname
	integer, intent(in) :: kind_of
	integer, parameter :: u = 42
	integer :: ib, is, k

	open(unit=u, file=fname, status='replace', action='write')

	!header and distinct bond lengths
	write(u,'(a)') '# File generated by extract_TB2J from '//infile
	write(u,'(a,f10.4,a)') '# Bond length cutoff = ', cutoff, ' A'
	if (home) then
		write(u,'(a)') '# Only ion pairs with R = (0, 0, 0)'
	else
		write(u,'(a)') '# Ion pairs with all R'
	end if
	write(u,'(i0)') nshell
	do is = 1, nshell
		write(u,'(f11.4)') shell(is)
	end do

	do is = 1, nshell
		do ib = 1, nbond
			if (bond_id(ib) /= is) cycle

		!line 1: ion_i ion_j bond_ID vector signs
		write(u,'(2x,a6,1x,a6,i6,3f10.4,3x,3i3)') &
			element_name(lab1(ib)), element_name(lab2(ib)), bond_id(ib), &
			vec(:,ib), (sign_of(vec(k,ib)), k=1,3)

		! next line(s): the value(s)
		select case (kind_of)
		case (1)	! J_iso
			write(u,'(2x,f12.6)') jiso(ib)
		case (2)	! J_ani
			do k = 1, 3
				write(u,'(2x,3f12.6)') jani(k,1,ib), jani(k,2,ib), jani(k,3,ib)
			end do
		case (3)	! DMI
			write(u,'(2x,3f12.6)') dmi(:,ib)
		case (4)	! Jprime
			write(u,'(2x,f12.6)') jprime(ib)
		end select
		end do
	end do

	close(u)
	end subroutine write_file

	end subroutine extract_TB2J
