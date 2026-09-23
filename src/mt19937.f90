	module mt19937
		use init, only: dp
		implicit none

		private
		integer, parameter :: i8 = 8
		integer, parameter :: nn = 624
		integer, parameter :: mm = 397
		integer(i8), parameter :: matrix_a   = int(z'9908b0df', i8)
		integer(i8), parameter :: upper_mask = int(z'80000000', i8)
		integer(i8), parameter :: lower_mask = int(z'7fffffff', i8)
		integer(i8), parameter :: full_mask  = int(z'ffffffff', i8)
		integer(i8), save :: mt(0:nn-1)
		integer, save :: mti = nn + 1
		public :: mt_init, mt_rand

	contains

	subroutine mt_init(seed)

		!  Initialize (or re-initialize) the generator state from an integer seed.
		integer, intent(in) :: seed
		integer :: i

		mt(0) = iand(int(seed, i8), full_mask)
		do i = 1, nn - 1
			mt(i) = iand(1812433253_i8 * ieor(mt(i-1), ishft(mt(i-1), -30)) + i, full_mask)
		end do
		mti = nn
	end subroutine mt_init

	function mt_rand() result(r)
		real(dp) :: r
		integer :: kk
		integer(i8) :: y
		integer(i8) :: mag01(0:1)

		mag01(0) = 0_i8
		mag01(1) = matrix_a

		if (mti >= nn) then

		! default seed if user forgot to call mt_init
		if (mti == nn + 1) call mt_init(5489)

		do kk = 0, nn - mm - 1
			y = ior(iand(mt(kk), upper_mask), iand(mt(kk+1), lower_mask))
			mt(kk) = ieor(ieor(mt(kk+mm), ishft(y, -1)), mag01(iand(y, 1_i8)))
		end do
		do kk = nn - mm, nn - 2
			y = ior(iand(mt(kk), upper_mask), iand(mt(kk+1), lower_mask))
			mt(kk) = ieor(ieor(mt(kk+mm-nn), ishft(y, -1)), mag01(iand(y, 1_i8)))
		end do
		y = ior(iand(mt(nn-1), upper_mask), iand(mt(0), lower_mask))
		mt(nn-1) = ieor(ieor(mt(mm-1), ishft(y, -1)), mag01(iand(y, 1_i8)))
		mti = 0
		end if

		y = mt(mti)
		mti = mti + 1

		! tempering
		y = ieor(y, ishft(y, -11))
		y = iand(ieor(y, iand(ishft(y, 7),  int(z'9d2c5680', i8))), full_mask)
		y = iand(ieor(y, iand(ishft(y, 15), int(z'efc60000', i8))), full_mask)
		y = ieor(y, ishft(y, -18))

		! divide by 2^32 -> [0,1)
		r = real(y, dp) / 4294967296.0_dp
	end function mt_rand

	end module mt19937
