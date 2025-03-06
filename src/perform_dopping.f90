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

	subroutine perform_dopping
	
	use init
	
	implicit none

	integer :: eof, IDs_to_interchange(2,nspecies), dopant, &
		ID_to_dope, dopant_ID, i, j, k, l, &
		available_sites_to_dope, total_sites_to_dope, ith_dope

	real(dp) :: dopants_spin, doping_ratio, dopants_spin_val(maxium_dopants)
	real(dp), allocatable :: s_dummy(:)

	character(len=80) :: text
	character(len=2) :: dope_at_species, dope_with_species, dopants_species(maxium_dopants)
	character(len=2), allocatable :: species_dummy(:)

	logical :: file_present

	dopants_species = ''
	dopants_spin_val = real(0.0, dp)

	inquire(file='dopants', exist=file_present)

	if(.not.file_present) then
	if (rank == 0) then
		write(6, *) "==> Required file 'dopants' is not present"
		write(6, *) "	 STOPPING now"
		write(6, *) ""
	end if
		stop
	end if

	open(unit=1, file='dopants', status='old', action='read')
	if (rank == 0) write(6, *) "==> Reading list of dopants"

	dopant = nspecies
	IDs_to_interchange = 0
	ndope = 0
	reading_dopants_list: do
	read(1, '(a)', end=10) text

	read(text, *) dope_at_species, dope_with_species, dopants_spin, doping_ratio

	if(any(dope_at_species.ne.species(1:))) then
		if (rank == 0) then
			write(6, *) "==> Magnetic ion at which doping will be performed"
			write(6, *) "    should be present in the 'structure.vasp' file."
			write(6, *) "	 STOPPING now"
			write(6, *) ""
		end if
		stop
	end if

	doping_ratio = abs(doping_ratio)

	if(doping_ratio.le.1) then

		ID_to_dope = findloc(species, value=dope_at_species, dim=1)
		ID_to_dope = ID_to_dope - 1
		IDs_to_interchange(1, ID_to_dope) = ID_to_dope
		if(any(dope_with_species.eq.species(1:))) then
			dopant_ID = findloc(species, value=dope_with_species, dim=1)
			dopant_ID = dopant_ID - 1
			IDs_to_interchange(2, ID_to_dope) = dopant_ID
		else
			dopant = dopant + 1	! Dopants ID
			ndope = ndope + 1	! Number of dopants
			dopants_species(dopant) = dope_with_species
			dopants_spin_val(dopant) = dopants_spin
			IDs_to_interchange(2, ID_to_dope) = dopant
		end if

		available_sites_to_dope = 0
		do l = 1, lattice_per_unit_cell
			do k = fromz, toz
				do j = fromy, toy
					do i = fromx, tox

		if(int(ion(4, i, j, k, l)).eq.IDs_to_interchange(1, ID_to_dope)) then
			available_sites_to_dope = available_sites_to_dope + 1
		end if
					end do
				end do
			end do
		end do

		total_sites_to_dope = nint(doping_ratio*available_sites_to_dope)
		if (rank == 0) write(6, "(' ==> Total sites for dopping: ',i6)") total_sites_to_dope
		ith_dope = 0
		dopping_sites: do
			call get_random_indices(i, j, k, l)
			if(int(ion(4, i, j, k, l)).eq.IDs_to_interchange(1, ID_to_dope)) then
				ion(4, i, j, k, l) = IDs_to_interchange(2, ID_to_dope)
				ith_dope = ith_dope + 1
			end if
			if(ith_dope.eq.total_sites_to_dope) exit dopping_sites

		end do dopping_sites

	else

	if (rank == 0) then
		write(6, *) "==> 'doping_ratio' cannot be greator than 1"
		write(6, *) "	  STOPPING now"
		write(6, *) ""
	end if
		stop
	end if

	end do reading_dopants_list
10	close(1)

	if(dopant.gt.nspecies) then

		allocate(s_dummy(nspecies), species_dummy(0:nspecies))
		species_dummy = ''
		s_dummy = s
		species_dummy = species

		deallocate(s, species, ionn)
		allocate(s(dopant), species(0:dopant), ionn(0:dopant))

		species(0) = species_dummy(0)
		do i = 1, dopant
			if(i.le.nspecies)then
				s(i) = s_dummy(i)
				species(i) = species_dummy(i)
			else
				s(i) = dopants_spin_val(i)
				species(i) = dopants_species(i)
			end if
		end do

		nspecies = dopant

		ionn = 0

		do l = 1, lattice_per_unit_cell
			do k = fromz, toz
				do j = fromy, toy
					do i = fromx, tox

		ionn(int(ion(4, i, j, k, l))) = ionn(int(ion(4, i, j, k, l))) + 1

					end do
				end do
			end do
		end do

		deallocate(s_dummy, species_dummy)
	end if


	end subroutine perform_dopping
