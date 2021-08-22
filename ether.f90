!				#################################
!				#	  ETHER			#
!				#	~ by Mukelian		#
!				#################################
																
!    Code Ether, based on Monte Carlo technique, can be used to study the 
!    thermodynamics of magnetic system for any crystal system.
!    Copyright (C) 2021 2021  Mukesh Kumar Sharma, Department of Physics,
!    Indian Institute of Technology Roorkee, Uttrakhand, India, PIN code 247667

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.

        program ether	!UPDATED Aug 22 2021 @ 02:54 PM
        implicit none
        integer, parameter :: dp = selected_real_kind(15,300)
        real(dp), allocatable :: a(:,:), ar(:,:)				!Atom
        integer, allocatable ::  nn(:,:,:,:), ionn(:)
        real(dp) :: ion_ref(5)							!trial ion
        real(dp), allocatable :: x(:), y(:), z(:)				!x,y,z co-ordinates of ions
        real(dp), allocatable :: s(:), sp(:,:), nbd_dis(:)
        real(dp), allocatable :: x_r(:), y_r(:), z_r(:)				!x,y,z co-ordinates of refrence ions
        real(dp), allocatable :: j_exc(:,:,:,:)
        integer, allocatable :: probe(:),t_probe(:), sc(:)			!to probe the special ions which will be used for only calculations
        real(dp), allocatable :: ion(:, :, :, :, :), ion_dummy(:, :, :, :, :)
	integer :: n_atoms, n_speci_incl,n_spc_excl, t_nbd_count, t_count
	integer :: atom_count, o_c, nbd_count, sc_i, s_count, total_ions
	real(dp) :: wt, abc(3,3), Si_dot_Sj(3), phi, s_trial_present(3)
	integer :: i, j, k, bulk_c, bond_c, sp_dim, n_species, no_of_nbd
	integer :: i_p, n_p, n_eq, accept_count, i_s, ii, jj, kk, ll
	integer :: stg
	real(dp) :: h_t, l_t, t_int !temperature :: High , Low, Interval
	integer :: sc_a, sc_b, sc_c ! supercell along a, b, c
	real(dp) :: nbd_range, dist_x, dist_y, dist_z, distance
	integer :: m, n, l, l_ion, sample
	real(dp) :: temp, beta, eta, sample_1, volume, h(3), g_factor
        real(dp) :: mag(3), rn, d_angle,  convert_to_rad, d_phi, la, lb, lc
	real(dp) :: kb = 8.6173303e-5_dp, pi = 3.141592653589793238462643383279_dp, mb
	character :: bc_x, bc_y, bc_z
        integer :: nscan, num, iscan, to_cal
	character*30 :: title, coordinate, filename
	character(len=2), allocatable :: species(:), speci_symb(:)
	logical :: opt_stg, field_b
	integer :: line
	real(dp) :: j_value(3), j_value_, anisotropy(3)

	character(len=20), dimension(50) :: m_head, e_head
        character*7 :: cc(7)
        
	real(dp) :: s_mag_avg, s_mag2_avg, err_mag_avg, err_mag2_avg
        real(dp) :: e_mag_avg,e_mag2_avg, mag_eng, trial_mag_eng
	real(dp) :: s_eng_avg, e_eng_avg, s_eng2_avg, eng
	real(dp) :: e_eng2_avg, err_eng_avg, err_eng2_avg
	real(dp) :: cv, s_cv, err_cv, e_cv
	real(dp) :: chi, s_chi, err_chi, e_chi
	real(dp) :: s_U_eng, U_eng, s_U_mag, U_mag
	real(dp) :: err_U_mag, err_U_eng, e_U_mag, e_U_eng
	real(dp) :: b_eng4_avg, mag4_avg
	real(dp) :: eng4_avg, eng_avg, eng2_avg, eng_tp
	real(dp) :: mag2_avg, mag_avg, mag_value
	
        real(dp) :: start, finish
        character(8)  :: date
        character(10) :: time
        character(5)  :: zone
        integer,dimension(8) :: values
        integer :: days, hrs, mins, secs, n_seed

	integer, dimension(:), allocatable :: seed
	real(dp) :: r
	
	call random_seed(size=n_seed)
	allocate(seed(n_seed))

        call cpu_time(start)
        call date_and_time(date,time,zone,values)
        call date_and_time(DATE=date,ZONE=zone)
        call date_and_time(TIME=time)
        call date_and_time(VALUES=values)

	call system ('rm -rf _spin _data *.dat')
	call system ('mkdir _spin _data')

	m_head(1) = '# 1 Temp.'
	m_head(2) =' 2 avg_mag'; m_head(3)='3 Chi'
	m_head(4) ='4 err_mag'; m_head(5) ='5 err_chi.'
	m_head(6) ='6 U_mag.'; m_head(7) ='7 err_U'

        e_head(1) = '# 1 Temp.'
	e_head(2) =' 2 avg_eng'; e_head(3)='3 Cv'
	e_head(4) ='4 err_eng'; e_head(5) ='5 err_Cv'
	e_head(6) ='6 U_eng'; e_head(7) ='7 err_U'

	cc(1) = ' # Atom '
	cc(2) = '    x';cc(3) = '  y';cc(4) = '  z'
	cc(5) = '  S_x';cc(6) = '  S_y';cc(7) = '  S_z'
		
	
	print*, "________________ Mukelian ________________"
        print*, ""
        print '(" PROGRAM STARTED on date ",i2,"-",i2,"-",i4)', values(3),values(2), values(1)
        print "(' at time ',i2,' hrs. ',i2,' min. ',i2,' sec. ')", values(5:7)
	print*,""
	print*,"Based on classical Heisenberg model; H = JSi*Sj"
	print*,""
	open(10001, file='structure.vasp', status='unknown')
	open(10003, file='initial_spin_conf.dat', status='unknown')	! initial structure file with intial spin
	open(unit=10004, file='magnetization.dat', status='replace'&
        ,action='write')
        open(unit=10005, file='energy.dat', status='replace'&
        ,action='write')
	open(10006, file='nbd.dat', status='unknown')
        open(unit=0, file='input', status='old', action='read')

	write(10004,10300) (m_head(i), i = 1, 7) 
	write(10005, 10300) (e_head(i), i = 1, 7)
	
	!______________ READING INPUT  _________________!
	read(0, *) n_p							! no. of pass
        read(0, *) n_eq							! equilibration steps
        read(0, *) h_t, l_t, t_int					! higher temp, lower temp., temp. interval
	read(0, *) n_species						! no. of species
	!# SPIN
	allocate(s(n_species))
	read(0, *) (s(i), i = 1, n_species)				! spin
	read(0, *) n_speci_incl						! no. of species to include
	allocate(speci_symb(n_speci_incl))
	read(0, *) (speci_symb(i), i = 1, n_speci_incl) 		! species lable to include in calculation
	read(0, *) sc_a, sc_b, sc_c					! supercell dimension
	read(0, *) opt_stg						! for staggered
	stg = 1								! default staggered value is 1
	read(0, *) bc_x, bc_y, bc_z					! boundary condition along x, y, z
	read(0, *) sample 						! sample per unit MC calculations)
	read(0, *) d_angle 						! least angle window, help in convergence
	read(0, *) to_cal 						! at this step interval all observables will be calculated
	read(0, *) field_b, h(1:3)					! Magnetic field logic, Mx, My, Mz (in T)
	read(0, *) g_factor						! g-factor
	close(0)							! closing of input file			@input

	!# EXCHANGE ENERGY
	open(10002, file='j_exchange_input', status='old', action='read')
	read(10002, *) no_of_nbd					! no. of diff. bond length w.r.t. any central ion
	allocate(nbd_dis(no_of_nbd))
	read(10002, *) nbd_dis(1:no_of_nbd)				! values of diff. bond lengths e.g., d1, d2, d3...
        allocate(j_exc(no_of_nbd, n_species, n_species,1:3))    	!j(:,:,xx-yy-zz)
	read(10002, *) line
	j_exc = 0.0_dp
	do ll = 1, line
		read(10002, *) i, j, k, j_value_, anisotropy(1:3)	! i = ith nbd; k, j is species comb.; corresponding j_value; it's anisotropy
			j_value = j_value_
			j_exc(i, j, k, 1:3) = j_value*anisotropy
			j_exc(i, k, j, 1:3) = j_exc(i, j, k, 1:3)	! making  sure for J(nbd1, nbd2), J(nbd2, nbd1) combination.
	end do

        	print*, "ALL values for 'J' should be provided in"
        	print*, "the terms of milli orders (meV), "
        	print*, "eg., for 1meV put only 1 not 0.001"
        	print*, ''
        	print*, '****************************************'
        	print*, "KINDLY CHECK THE 'j_exchange_input' FILE"
        	print*, '(IGNOR, if it is already DONE!)'
        	print*, '****************************************'
        	print*, ''

	j_exc = j_exc*(0.001_dp)/kb					! converting 'j_exc' into temp.

	close(10002) !closing of j_exchange file

        mb = 5.7883818060e-5_dp/kb					! converting bhor magneton in K/T

	!LEAST ANGLE FOR VARYING THE PHI ANGLE FOR CONVERGENCE
	convert_to_rad = pi/180.0_dp
	d_phi = convert_to_rad*d_angle
	
	! READING STRUCTURE FILE   ###############################################
	!#########################################################################
	read(10001,*) title
	read(10001,*) wt
	lattice_constant : do i = 1,3
	read(10001,*) (abc(i,j), j= 1,3)
	enddo lattice_constant
	la = sqrt(dot_product(abc(1, 1:3), abc(1, 1:3)))	! lattice parameter a
	lb = sqrt(dot_product(abc(2, 1:3), abc(2, 1:3)))	! lattice parameter b
	lc = sqrt(dot_product(abc(3, 1:3), abc(3, 1:3)))	! lattice parameter c
	allocate(ionn(0:n_species), species(n_species))
	read(10001,*) (species(i), i =1,n_species)
	ionn = 0
	read(10001,*) (ionn(i), i =1 ,n_species)
	n_atoms = sum(ionn)

	allocate(x(n_atoms), y(n_atoms), z(n_atoms))

	x = 0; y = 0; z = 0
	read(10001,*) coordinate
	if(coordinate.eq.'Direct')then
		print*, 'Choose POSCAR in cartesian co-ordinate only'
		print*, 'STOPPING'
		stop
	endif

	do i = 1,n_atoms
		read(10001,*) x(i), y(i), z(i) !cartesian co-ordinates
	end do
	close(10001)
	print*, 'Reading structure is completed!'
	!#########################################################################

	!SAVING DETAILS OF ATOMS  ################################################
	!#########################################################################
	allocate(a(n_atoms,0:8))
	a = 0
	atomic_details: do k = 1, n_species
					do i = int(sum(ionn(0:k-1)))+1, int(sum(ionn(0:k)))

						a(i,0) = i								! atom no.
						call random_uniform(-1.0_dp, +1.0_dp, rn)
						a(i, 3) = rn								! S(z)
						call random_uniform(0.0_dp, 2.0_dp*pi, rn)
        					a(i, 1) = sqrt(1.0_dp-(a(i, 3)**2.0_dp))*cos(rn)	! S(x)
						a(i, 2) = sqrt(1.0_dp-(a(i, 3)**2.0_dp))*sin(rn)	! S(y)
						a(i, 4) = k								! tag for element
						a(i, 5) = rn								! phi value
						a(i, 6) = x(i)							! x
						a(i, 7) = y(i)							! y
						a(i, 8) = z(i)							! z
						a(i, 1:3) = (a(i, 1:3)/&						! normalization
						sqrt(dot_product(a(i, 1:3),a(i, 1:3))))*s(k)      
					enddo
				enddo atomic_details

        print'(3I4," :: (Lx, Ly, Lz) <==> super cell size")',sc_a, sc_b, sc_c
        print*,''
	print*, 'NOTE: Working on MCS per spin'
        print*,''
        print'(3A4," :: (along x,y,z-axis boundary condition)")',bc_x, bc_y, bc_z
        print*,"                'o' => open; 'c' => closed"
        print*,''
	print*,'==> spins are set in random configuration'
        print'(" ==> +/-",f6.2," deg. angle has been chosen for the convergence")', d_angle
	if(opt_stg)then
		print*,' ==> Staggered magnetization is SELECTED!'
	end if
	if(field_b)then
		print'(" ==> Magnetic field ::", 3f10.5, " (Mx, My, Mz) is SELECTED!")', h
	end if
	print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	print*, ''
	!#########################################################################
	!INCLUSION OF SPECIES FOR CALC.  #########################################
	!#########################################################################
	s_count = 0
	do i = 1, n_speci_incl
		do j =1 , n_species
			if (speci_symb(i).eq.species(j)) then
				s_count = s_count + ionn(j)
			end if
		end do
	end do

	allocate(ar(s_count, 0:8))
	
	atom_count = 0
	inclusion : do i = 1, n_speci_incl
		do j =1 , n_species
			if (speci_symb(i).eq.species(j)) then
				do k = sum(ionn(0:j-1))+1, sum(ionn(0:j))
					atom_count = atom_count +1	! Counting included ions
					ar(atom_count,:) = a(k,:)	! 'k' will be used for identifying
					 				! these element in further last moment 
				enddo
			endif
		enddo
	enddo inclusion
	!#########################################################################

	volume = dble(sc_a*sc_b*sc_c*atom_count)
        
	!EXPANDING INTO SUPERCELL  ##############################################
	!#########################################################################
	allocate(ion(sc_a+2, sc_b+2, sc_c+2, atom_count, 0:8), ion_dummy(sc_a+2, sc_b+2, sc_c+2, atom_count, 0:8))
	s_count = 0
	do k = 1, sc_c + 2
		do j = 1, sc_b + 2
			do i = 1, sc_a + 2
				do l = 1, atom_count
					s_count = s_count +1		! Counting total no. of ions in supercell
					ion(i, j, k, l, 0) = s_count	! ID
					ion(i, j, k, l, 1:5) = ar(l, 1:5)	! Sx, Sy, Sz, species no. & phi
					!ion(i, j, k, l, 6) = ar(l, 6) + la*(i -1)	! x co-ordinate
					!ion(i, j, k, l, 7) = ar(l, 7) + lb*(j -1) 	! y co-ordinate
					!ion(i, j, k, l, 8) = ar(l, 8) + lc*(k -1)	! z co-ordinate
					ion(i, j, k, l, 6:8) = abc(1, 1:3)*(i -1) &	! x co-ordinate
								+ abc(2, 1:3)*(j -1) + abc(3, 1:3)*(k -1) &
								+ ar(l, 6:8)
				end do
			end do
		end do
	end do
	total_ions = s_count
	!#########################################################################

	!BOUNDARY CONDITION  ###################################################
	!#########################################################################	

	! BOUNDAR Z
	boundary_zo : if(bc_z.eq.'c')then
					ion(:, :, 1, :, 1:5) = ion(:, :, sc_c +1, :, 1:5)
					ion(:, :, sc_c +2, :, 1:5) = ion(:, :, 2, :, 1:5)												
				else if (bc_z.eq.'o')then				
					ion(:, :, 1, :, 4) = 0
					ion(:, :, sc_c +2, :, 4) = 0			
				else	
					print*, 'ERROR:: Check the boundary condition input on "Z"'
					stop
																				
				end if boundary_zo
								
	! BOUNDAR Y
	boundary_yo : if(bc_y.eq.'c')then		
					ion(:, 1, :,  :, 1:5) = ion(:, sc_b +1, :, :, 1:5)
					ion(:, sc_b +2,  :, :, 1:5) = ion(:, 2, :,  :, 1:5)					
				elseif (bc_y.eq.'o')then
					ion(:, 1, :,  :, 4) = 0
					ion(:, sc_b +2,  :, :, 4) = 0				
				else	
					print*, 'ERROR:: Check the boundary condition input on "Y"'
					stop
																	 											 
				end if boundary_yo
								
	! BOUNDAR X
	boundary_xo : if(bc_x.eq.'c')then			
					ion(1, :, :,  :, 1:5) = ion(sc_a +1, :, :, :, 1:5)
					ion(sc_a +2, :,  :, :, 1:5) = ion( 2, :, :,  :, 1:5)					
				elseif (bc_x.eq.'o')then
					ion(1, :, :,  :, 4) = 0
					ion(sc_a +2, :,  :, :, 4) = 0				
				else	
					print*, 'ERROR:: Check the boundary condition input on "X"'
					stop
											
				end if boundary_xo

	!NEAREST NEIGHBOUR   ####################################################
	!#########################################################################
	allocate(nn(no_of_nbd, total_ions, 0:20, 0:4))
		nn = 0
		central_atom : do k = 1, sc_c + 2
				do j = 1, sc_b + 2
					do i = 1, sc_a + 2
							do l = 1, atom_count
								nbd : do m = 1, no_of_nbd	! no. of distinct bond
								nbd_count = 0.0
							
				do kk = 1, sc_c + 2
					do jj = 1, sc_b + 2
						do ii = 1, sc_a + 2
							do ll = 1, atom_count

							                if(ion(i, j, k, l, 0).ne.ion(ii, jj, kk, ll, 0))then                 
											distance = sqrt(dot_product(&
											ion(i, j, k, l, 6:8) - ion(ii, jj, kk, ll, 6:8), &
											ion(i, j, k, l, 6:8) - ion(ii, jj, kk, ll, 6:8)  &
											))
										
											if(abs(nbd_dis(m) - (distance)).le.0.001)then
												nbd_count = nbd_count +1							
												nn(m, int(ion(i, j, k, l, 0)), nbd_count, 0) = ion(ii, jj, kk, ll, 0)	! Storing nbd's ID
												nn(m, int(ion(i, j, k, l, 0)), nbd_count, 1) = ii			! Storing nbd's indices
												nn(m, int(ion(i, j, k, l, 0)), nbd_count, 2) = jj			! Storing nbd's indices
												nn(m, int(ion(i, j, k, l, 0)), nbd_count, 3) = kk			! Storing nbd's indices
												nn(m, int(ion(i, j, k, l, 0)), nbd_count, 4) = ll			! Storing nbd's atom count
                                                                                        end if        
                                                                        end if

							end do
						end do
					end do
				end do
									nn(m, int(ion(i, j, k, l, 0)), 0, 0) = nbd_count	! Storing total no. of similar
																! nbds on particular ion for ith distinct bond
								end do nbd								
							end do
						end do
					end do
				end do central_atom
	!#########################################################################

	! WRITING DATAS FOR NBD  #################################################
	!#########################################################################
	nbd_inf : do k = 1, sc_c + 2
				do j = 1, sc_b + 2
					do i = 1, sc_a + 2
                                                do l = 1, atom_count
							write(10006,10) int(ion(i, j, k, l, 0))
							write(10006,*)"~~~~~~~~~~~~~~" 
							nbd_ : do m = 1, no_of_nbd
								write(10006, 11) m, &
								nn(m, int(ion(i, j, k, l, 0)), &
								1:nn(m, int(ion(i, j, k, l, 0)), 0, 0), &
								 0)
10      format(" ION no. ",i5)
11      format(" For bond length no. ",i2," nbds are ==> "20i5)
							end do nbd_
                                                        write(10006,*)''
                                                end do
					end do
				end do
		end do nbd_inf
	!#########################################################################
				
	! WRITING INITIAL STRUCTURE FILE   #########################################
	!#########################################################################
	write(10003,10030) (cc(j),j=1, 7)	! for xcrysden file formate 
	write(10003,*) 'ATOM'
			do k = 1, sc_c + 2
				do j = 1, sc_b + 2
					do i = 1, sc_a + 2
						do l = 1, atom_count
							write(10003,10021) species(int(ion(i, j, k, l, 4))), ion(i, j, k, l, 6:8), ion(i, j, k, l, 1:3)
						end do
					end do
				end do
			end do
	close(10003)
	!#########################################################################

        !SCAN LOOP ##############################################################
       	!#########################################################################
       	num = 20000
        nscan = int((h_t - l_t)/(t_int)) + 1
        scan_loop : do iscan = 1, nscan	
        temp = h_t - t_int*(iscan-1)

        !FORMATION OF FILES @ temp. K  ############################################
       	!#########################################################################
        num = num + 1
        write(filename,100) temp
100   	format('sp', f7.2, 'K.dat')
        open(file = filename, unit = num)
       	!#########################################################################
        
        write(num,10030)(cc(j),j=1, 7)	! for xcrysden file formate 
	write(num,*) 'ATOM'
   
        beta = 1.0_dp/temp

        !SAMPLING  ###############################################################
       	!#########################################################################
       	
	print "(' Started for temp.', f10.5, ' K')",temp
	print "(' ``````````````````````````````')"
	print*,''
        
	s_eng_avg = 0.0_dp; e_eng_avg = 0.0_dp
	s_eng2_avg = 0.0_dp; e_eng2_avg = 0.0_dp
	s_U_eng = 0.0_dp; e_U_eng = 0.0_dp
	s_cv = 0.0_dp; e_cv = 0.0_dp

	s_mag_avg = 0.0_dp; e_mag_avg = 0.0_dp
	s_mag2_avg = 0.0_dp; e_mag2_avg = 0.0_dp
	s_U_mag = 0.0_dp; e_U_mag = 0.0_dp
	s_chi = 0.0_dp; e_chi = 0.0_dp

       	ion_ref = 0.0_dp; trial_mag_eng = 0.0; mag_eng = 0.0_dp
	
        sampling : do i_s = 1, sample
        
        eng_avg = 0.0_dp; eng2_avg = 0.0_dp; eng4_avg = 0.0_dp
        mag_avg = 0.0_dp; mag2_avg = 0.0_dp; mag4_avg = 0.0_dp
        o_c = 0.0_dp
        
        seed = i_s
        call random_seed(put=seed)	! feeding new seed for new series of random numbers
        
        accept_count = 0.0
        if (iscan.ne.1)then
        	ion = ion_dummy	!set the spins to optimized spins done at earlier temp. (t-dt) K.
	end if

        !MONTE CARLO  ###########################################################
       	!#########################################################################

        Monte_Carlo : do i_p = 1, n_p
        
	!MCS PROCESS ###########################################################
       	!#########################################################################

        MCS_calculation : do k = 2, sc_c +1
					do j = 2, sc_b +1
						do i = 2, sc_a +1
							do l = 1, atom_count
								ion_ref(4) = ion(i, j, k, l, 4)						! storing species no.
								ion_ref(5) = ion(i, j, k, l, 5)						! old phi
								call random_uniform(-1.0_dp, +1.0_dp, rn)
								ion_ref(3) = rn							! cos(theta) : S(z)
								call random_uniform(-1.0_dp, +1.0_dp, rn)
								phi = ion_ref(5) + (d_phi*rn)						!older_phi+ delta*rn
        							ion_ref(1) = sqrt(1.0_dp-(ion_ref(3)**2.0_dp))*cos(phi)	! S(x)
								ion_ref(2) = sqrt(1.0_dp-(ion_ref(3)**2.0_dp))*sin(phi)	! S(y)

								ion_ref(1:3) = (ion_ref(1:3)/sqrt(dot_product(ion_ref(1:3),ion_ref(1:3))))*s(int(ion_ref(4)))
								ion_ref (5) = phi							! storing phi value
								eng_tp = 0.0_dp								! energy trial present
								
								s_trial_present = ion_ref(1:3) - ion(i, j, k, l, 1:3)			! S_trial - S_present
	
									do m = 1, no_of_nbd						! for distinct bond
										do n = 1, nn(m, int(ion(i, j, k, l, 0)), 0, 0)		! no. of similar nbd for ith distinct bond

											if(ion(nn(m, int(ion(i, j, k, l, 0)), n, 1), &
											nn(m, int(ion(i, j, k, l, 0)), n, 2), &
											nn(m, int(ion(i, j, k, l, 0)), n, 3), &
											nn(m, int(ion(i, j, k, l, 0)), n, 4), 4).eq.0) then
											go to 1
											end if

											Si_dot_Sj = s_trial_present &			!Si*Sj
											* ion(nn(m, int(ion(i, j, k, l, 0)), n, 1), &
											nn(m, int(ion(i, j, k, l, 0)), n, 2), &
											nn(m, int(ion(i, j, k, l, 0)), n, 3), &
											nn(m, int(ion(i, j, k, l, 0)), n, 4), 1:3)

											eng_tp = eng_tp + &
												dot_product(j_exc(m, int(ion(i, j, k, l, 4)), &
												int(ion(nn(m, int(ion(i, j, k, l, 0)), n, 1), &
												nn(m, int(ion(i, j, k, l, 0)), n, 2), &
												nn(m, int(ion(i, j, k, l, 0)), n, 3), &
												nn(m, int(ion(i, j, k, l, 0)), n, 4), 4)), &
												  :), &
												Si_dot_Sj)

1												continue
										end do
									end do

								magnetic_tp : if(field_b)then !magnetic part
									trial_mag_eng = -(g_factor*mb*dot_product(s_trial_present, h))	! energy due to magnetic field
									eng_tp = eng_tp + trial_mag_eng
								end if magnetic_tp

								!METROPOLIS ALGO #################################
								!###################################################
								call random_uniform(0.0_dp,1.0_dp,rn)
								eta = rn ! random number 0-1
								metropolis : if (exp(-beta*eng_tp) .gt. eta) then	! exp(-deltaE(trial - present)*beta) > eta
											ion(i, j, k, l, 1:5) = ion_ref(1:5)
											accept_count = accept_count +1.0
								
									! BOUNDAR Z
									boundary_z : if(bc_z.eq.'c')then
													if (k.eq.(sc_c +1)) then
														ion(i, j, 1, l, 1:5) = ion(i, j, k, l, 1:5)
													elseif(k.eq.2)then
														ion(i, j, sc_c +2, l, 1:5) = ion(i, j, k, l, 1:5)
													end if
												end if boundary_z
								
									! BOUNDAR Y
									boundary_y : if(bc_y.eq.'c')then
													if (j.eq.(sc_b +1)) then
														ion(i, 1, k, l, 1:5) = ion(i, j, k, l, 1:5)
													elseif(j.eq.2)then
														ion(i, sc_b +2, k, l, 1:5) = ion(i, j, k, l, 1:5)
													end if
												end if boundary_y
								
									! BOUNDAR X
									boundary_x : if(bc_x.eq.'c')then
													if (i.eq.(sc_a +1)) then
														ion(1, j, k, l, 1:5) = ion(i, j, k, l, 1:5)
													elseif(i.eq.2)then
														ion(sc_a +2, j, k, l, 1:5) = ion(i, j, k, l, 1:5)
													end if
												end if boundary_x
										end if metropolis
								!###################################################
								
							end do
						end do
					end do
				end do MCS_calculation

        equilibration : if((i_p.gt.n_eq).and.(mod(real(i_p), real(to_cal)).eq.0.0))then

        o_c = o_c + 1.0_dp
	mag = 0.0_dp ; eng = 0.0_dp

        calculation : do k = 2, sc_c +1
				do j = 2, sc_b +1
					do i = 2, sc_a +1
						do l = 1, atom_count
							if(opt_stg.eqv..TRUE.) stg = (-1)**(i+j+k+l)
							
								mag = mag + stg*ion(i, j, k, l, 1:3)	! magnetization vector
								
								do m = 1, no_of_nbd						! for distinct bond
									do n = 1, int(nn(m, int(ion(i, j, k, l, 0)), 0, 0))	! no. of similar nbd for ith distinct bond
									
										if(ion(nn(m, int(ion(i, j, k, l, 0)), n, 1), &
										nn(m, int(ion(i, j, k, l, 0)), n, 2), &
										nn(m, int(ion(i, j, k, l, 0)), n, 3), &
										nn(m, int(ion(i, j, k, l, 0)), n, 4), &
										 4).eq.0) then
										 go to 2
										 end if

										Si_dot_Sj = ion(i, j, k, l, 1:3) &		!Si*Sj
										* ion(nn(m, int(ion(i, j, k, l, 0)), n, 1), &
										nn(m, int(ion(i, j, k, l, 0)), n, 2), &
										nn(m, int(ion(i, j, k, l, 0)), n, 3), &
										nn(m, int(ion(i, j, k, l, 0)), n, 4), 1:3)

										eng = eng + &
											dot_product(j_exc(m, int(ion(i, j, k, l, 4)), &
											int(ion(nn(m, int(ion(i, j, k, l, 0)), n, 1), &
											nn(m, int(ion(i, j, k, l, 0)), n, 2), &
											nn(m, int(ion(i, j, k, l, 0)), n, 3), &
											nn(m, int(ion(i, j, k, l, 0)), n, 4), 4)), &
											  :), &
											Si_dot_Sj)
2											continue
									end do
									end do
								
						end do
					end do
				end do
			end do calculation

			magnetic_ : if(field_b)then !magnetic part
				mag_eng = -(g_factor*mb*dot_product(mag, h))	! energy due to magnetic field
				eng = eng + mag_eng
			end if magnetic_

	mag_value =  sqrt(dot_product(mag, mag))
        mag_avg = mag_avg + mag_value
        mag2_avg = mag2_avg + mag_value**2
        mag4_avg = mag4_avg + mag_value**4
        
        eng_avg = eng_avg + eng
        eng2_avg = eng2_avg + eng**2
        eng4_avg = eng4_avg + eng**4 
       
        end if equilibration

	end do Monte_Carlo
	!#########################################################################

	!AVG. MAGNETIZATION
	!################################################################
	mag_avg = mag_avg/(volume*o_c); mag2_avg = mag2_avg/((volume**2)*o_c)
	mag4_avg = mag4_avg/((volume**4)*o_c)
	s_mag_avg = mag_avg + s_mag_avg!; s_mag2_avg = mag2_avg + s_mag2_avg
	U_mag = 1.0_dp - ((1.0_dp/3.0_dp)*(mag4_avg/(mag2_avg**2)))
	s_U_mag = U_mag + s_U_mag
	e_U_mag = U_mag**2 + e_U_mag
	e_mag_avg = (mag_avg**2) + e_mag_avg!; e_mag2_avg = (mag2_avg**2) + e_mag2_avg
	chi = beta*(mag2_avg - (mag_avg**2))*volume
	s_chi = chi + s_chi; e_chi = (chi**2) + e_chi
	!################################################################

	!AVG. ENERGY
	!################################################################
	eng_avg = eng_avg/(2*volume*o_c); eng2_avg = eng2_avg/(((2*volume)**2)*o_c)
	eng4_avg = eng4_avg/(((2*volume)**4)*o_c)
	s_eng_avg = eng_avg + s_eng_avg!; s_eng2_avg = eng2_avg + s_eng2_avg
	U_eng = 1.0_dp - ((1.0_dp/3.0_dp)*(eng4_avg/(eng2_avg**2)))
	s_U_eng = U_eng + s_U_eng
	e_U_eng = e_U_eng + U_eng**2
	e_eng_avg = (eng_avg**2) + e_eng_avg!; e_eng2_avg = (eng2_avg**2) + e_eng2_avg
	cv = (beta*beta)*(eng2_avg - (eng_avg**2))*volume
	s_cv = cv + s_cv; e_cv = (cv**2) + e_cv
	!################################################################
	
	!################################################################
        print'(" completed for sample ID ",i2," with acceptance probability",&
        f6.2," %")',i_s,accept_count/(volume*n_p)*100

	end do sampling
	
	ion_dummy = ion	!storing all optimized spins into sample_array
	print'(" DONE!")'
	print*,''
	!################################################################

	sample = sample*1.0_dp
	sample_1 = 1.0_dp/(sample - 1.0_dp)
	
	!MAGNETIZATION PER SAMPLE
	!################################################################
	s_mag_avg = s_mag_avg/sample!; s_mag2_avg = s_mag2_avg/sample
	e_mag_avg = e_mag_avg/sample!; e_mag2_avg = e_mag2_avg/sample
	s_U_mag = s_U_mag/sample; e_U_mag = e_U_mag/sample
	err_U_mag = sqrt(sample_1*(e_U_mag - (s_U_mag**2)))
	s_chi = s_chi/sample; e_chi = e_chi/sample
	err_mag_avg = sqrt(sample_1*(e_mag_avg - (s_mag_avg**2)))
	!err_mag2_avg = sqrt(sample_1*(e_mag2_avg - (s_mag2_avg**2)))
	err_chi = sqrt(sample_1*(e_chi - (s_chi**2)))
	!################################################################
	
	!ENERGY PER SAMPLE
	!################################################################
	s_eng_avg = s_eng_avg/sample!; s_eng2_avg = s_eng2_avg/sample	
	e_eng_avg = e_eng_avg/sample!; e_eng2_avg = e_eng2_avg/sample
	s_U_eng = s_U_eng/sample; e_U_eng = e_U_eng/sample
	err_U_eng = sqrt(sample_1*(e_U_eng - (s_U_eng**2)))
	s_cv = s_cv/sample; e_cv = e_cv/sample
	err_eng_avg = sqrt(sample_1*(e_eng_avg - (s_eng_avg**2)))
	!err_eng2_avg = sqrt(sample_1*(e_eng2_avg - (s_eng2_avg**2)))	
	err_cv = sqrt(sample_1*(e_cv - (s_cv**2)))
	!################################################################
	
	!WRITTING OUTPUT FILES  ##################################################
	!#########################################################################
	
	!SPINS
			do k = 2, sc_c + 1
				do j = 2, sc_b + 1
					do i = 2, sc_a + 1
						do l = 1, atom_count
							write(num, 10021) species(int(ion(i, j, k, l, 4))), ion(i, j, k, l, 6:8), ion(i, j, k, l, 1:3)
						end do
					end do
				end do
			end do

	close(num)
	!#########################################################################

	!ENERGY
        write(10004,10200)temp, s_mag_avg, &
	s_chi, err_mag_avg, &
	err_chi, s_U_mag, err_U_mag
        
        !MAGNETIC
        write(10005,10200)temp, s_eng_avg, &
	s_cv, err_eng_avg, &
	err_cv, s_U_eng, err_U_eng

        end do scan_loop
        
        call cpu_time(finish)
        days = int((finish-start)/(24*3600.0))
        hrs = int(mod((finish-start), 24*3600.0)/3600.0)
        mins = int(mod(mod((finish-start), 24*3600.0), 3600.0)/60.0)
        secs = int(mod(mod(mod((finish-start), 24*3600.0), 3600.0), 60.0))

        print*,''
        print '(" Total calculation time = ",i3," days",i3," Hr",i3,&
        " min.",i3," sec.")',&
        days, hrs, mins, secs

        call system ('bash reduce.sh')
        call system ('mv sp* _spin')
        call system ('mv *.dat _data')
        call system ('cp graph* _data')

        print*,''
        print*,"May the force be with you"
        print*,"~ Mukesh Kumar Sharma"
        print*,"e-mail@ msharma1@ph.iitr.ac.in"
        print*,''
        close(10004)
        close(10005)

        call date_and_time(date,time,zone,values)
        call date_and_time(DATE=date,ZONE=zone)
        call date_and_time(TIME=time)
        call date_and_time(VALUES=values)
        print*,''
        print '(" PROGRAM ENDED on date ",i2,"-",i2,"-",i4)',values(3),values(2), values(1)
        print "(' at time ',i2,' hrs. ',i2,' min. ',i2,' sec. ')", values(5:7)
        print*,''

10030   format(12A7)
10200   format(50E21.9)
10021   format(A5,7f13.7)
10300   format(35A21)

        contains

        !_______RANDOM NUMBER GENERATING FUNCTION______!

        subroutine random_uniform(a,b,rn)
           implicit none
           integer, parameter :: dp = selected_real_kind(15,300)
           real(dp),intent(in) :: a,b
           real(dp),intent(out) :: rn
           real(dp) :: u
           call random_number(u)
           rn = (b-a)*u + a
        end subroutine random_uniform
        
        end program ether
