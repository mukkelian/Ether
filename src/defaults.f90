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

        subroutine defaults
        use init

        implicit none

        model = "Ising"
        Ising = .TRUE.
        tmcs = 5000
        tmcs_eq = 3000
        bc = 'c'
        ht = 100
        lt = 10
        tint = 5
        sc = (/2, 2, 2/)
        repeat = 10
        angle = .FALSE.
        to_cal = 20
        Zeeman = .FALSE.
        h = 0
        g_factor = 2
        SIA = .FALSE.
        para = .FALSE.
        J_para = 1.0
        ovrr = .FALSE.
        seed = 1992
        dope = .FALSE.

        end subroutine defaults
