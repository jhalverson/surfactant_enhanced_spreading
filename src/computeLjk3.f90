	subroutine computeLjk3(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ip(3, 8), kk, jj, ii, j, k, l, n, m, ib
	integer*4 icoeff, icoeff2, igetIndex, iIkm2
	real*8 rxbb, rybb, rzbb, rctr3x, rctr3y, rctr3z
	real*8 rrho, ralpha, rbeta, rAnmjk, rAnm, rtmp
	complex*16 ctmp, Ylm

	ip(1, 1) = 3
	ip(2, 1) = 3
	ip(3, 1) = 3

	ip(1, 2) = 4
	ip(2, 2) = 3
	ip(3, 2) = 3

	ip(1, 3) = 3
	ip(2, 3) = 4
	ip(3, 3) = 3

	ip(1, 4) = 4
	ip(2, 4) = 4
	ip(3, 4) = 3

	ip(1, 5) = 3
	ip(2, 5) = 3
	ip(3, 5) = 4

	ip(1, 6) = 4
	ip(2, 6) = 3
	ip(3, 6) = 4

	ip(1, 7) = 3
	ip(2, 7) = 4
	ip(3, 7) = 4

	ip(1, 8) = 4
	ip(2, 8) = 4
	ip(3, 8) = 4

	do ib = 1, 8

		rctr3x = rcen3x(ib)
		rctr3y = rcen3y(ib)
		rctr3z = rcen3z(ib)

		do kk = 1, 6
		do jj = 1, 6
		do ii = 1, 6

		if(abs(ii - ip(1, ib)).gt.1.or.abs(jj - ip(2, ib)).gt.1.or.abs(kk - ip(3, ib)).gt.1) then

			l = ii + (jj - 1) * 6 + (kk - 1) * 6 * 6
			rxbb = rcen3xLjk3(l) - rctr3x
			rybb = rcen3yLjk3(l) - rctr3y
			rzbb = rcen3zLjk3(l) - rctr3z

			rrho = sqrt(rxbb**2 + rybb**2 + rzbb**2)
			ralpha = acos(rzbb / rrho)
			rbeta = atan2(rybb, rxbb)

			do j = 0, p

				do k = 0, j

					icoeff = igetIndex(j, k)
					rAnmjk = rAnm(j, k)

					do n = 0, p

						do m = -n, n

						ctmp = cO3(igetIndex(n, m), l) * Ylm(j + n, m - k, ralpha, rbeta)
						rtmp = iIkm2(k, m) * rAnmjk * rAnm(n, m) / ((-1)**n * rAnm(j + n, m - k) * rrho**(j + n + 1))
						cL3(icoeff, ib) = cL3(icoeff, ib) + rtmp * ctmp

						enddo

					enddo

					if(k.gt.0) then

						icoeff2 = igetIndex(j, -k)
						cL3(icoeff2, ib) = conjg(cL3(icoeff, ib))

					endif

				enddo

			enddo

		endif

		enddo
		enddo
		enddo

	enddo

	return
	end subroutine
