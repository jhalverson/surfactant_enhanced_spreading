	subroutine computeLjk4(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 key4(8), key4i3, ib4(8), ip(3, 8), jp(3, 8)
	integer*4 i3, i4, ib4i4, ib512, ib64to512
	integer*4 ip1i3, ip2i3, ip3i3
	integer*4 ilowx, ilowy, ilowz, ihghx, ihghy, ihghz
	integer*4 kk, jj, ii, j, k, l, n, m
	integer*4 icoeff, icoeff2, igetIndex, iIkm2
        real*8 rxbb, rybb, rzbb
        real*8 rrho, ralpha, rbeta, rAnmjk, rAnm, rtmp
        complex*16 ctmp, Ylm

	key4(1) = 1
	key4(2) = 3
	key4(3) = 9
	key4(4) = 11
	key4(5) = 33
	key4(6) = 35
	key4(7) = 41
	key4(8) = 43

	ip(1, 1) = 3
	ip(2, 1) = 3
	ip(3, 1) = 3

	ip(1, 2) = 5
	ip(2, 2) = 3
	ip(3, 2) = 3

	ip(1, 3) = 3
	ip(2, 3) = 5
	ip(3, 3) = 3

	ip(1, 4) = 5
	ip(2, 4) = 5
	ip(3, 4) = 3

	ip(1, 5) = 3
	ip(2, 5) = 3
	ip(3, 5) = 5

	ip(1, 6) = 5
	ip(2, 6) = 3
	ip(3, 6) = 5

	ip(1, 7) = 3
	ip(2, 7) = 5
	ip(3, 7) = 5

	ip(1, 8) = 5
	ip(2, 8) = 5
	ip(3, 8) = 5

	do i3 = 1, 8

		key4i3 = key4(i3)

		ib4(1) = key4i3
		ib4(2) = key4i3 + 1
		ib4(3) = key4i3 + 4
		ib4(4) = key4i3 + 5
		ib4(5) = key4i3 + 16
		ib4(6) = key4i3 + 17
		ib4(7) = key4i3 + 20
		ib4(8) = key4i3 + 21

		ip1i3 = ip(1, i3)
		ip2i3 = ip(2, i3)
		ip3i3 = ip(3, i3)

		ilowx = ip1i3 - 2
		ihghx = ip1i3 + 3
		ilowy = ip2i3 - 2
		ihghy = ip2i3 + 3
		ilowz = ip3i3 - 2
		ihghz = ip3i3 + 3

		jp(1, 1) = ip1i3
		jp(2, 1) = ip2i3
		jp(3, 1) = ip3i3

		jp(1, 2) = ip1i3 + 1
		jp(2, 2) = ip2i3
		jp(3, 2) = ip3i3

		jp(1, 3) = ip1i3
		jp(2, 3) = ip2i3 + 1
		jp(3, 3) = ip3i3

		jp(1, 4) = ip1i3 + 1
		jp(2, 4) = ip2i3 + 1
		jp(3, 4) = ip3i3

		jp(1, 5) = ip1i3
		jp(2, 5) = ip2i3
		jp(3, 5) = ip3i3 + 1

		jp(1, 6) = ip1i3 + 1
		jp(2, 6) = ip2i3
		jp(3, 6) = ip3i3 + 1

		jp(1, 7) = ip1i3
		jp(2, 7) = ip2i3 + 1
		jp(3, 7) = ip3i3 + 1

		jp(1, 8) = ip1i3 + 1
		jp(2, 8) = ip2i3 + 1
		jp(3, 8) = ip3i3 + 1

		do i4 = 1, 8

			ib4i4 = ib4(i4)

			do kk = ilowz, ihghz
			do jj = ilowy, ihghy
			do ii = ilowx, ihghx

			if(abs(ii - jp(1, i4)).gt.1.or.abs(jj - jp(2, i4)).gt.1.or.abs(kk - jp(3, i4)).gt.1) then

				l = ii + (jj - 1) * 8 + (kk - 1) * 8 * 8

				rxbb = rcen4xLjk4(l) - rcen4x(ib4i4)
				rybb = rcen4yLjk4(l) - rcen4y(ib4i4)
				rzbb = rcen4zLjk4(l) - rcen4z(ib4i4)

				rrho = sqrt(rxbb**2 + rybb**2 + rzbb**2)
				ralpha = acos(rzbb / rrho)
				rbeta = atan2(rybb, rxbb)

				do j = 0, p

					do k = 0, j

						icoeff = igetIndex(j, k)
						rAnmjk = rAnm(j, k)

						do n = 0, p

							do m = -n, n

						ctmp = cO4(igetIndex(n, m), l) * Ylm(j + n, m - k, ralpha, rbeta)
						rtmp = iIkm2(k, m) * rAnmjk * rAnm(n, m) / ((-1)**n * rAnm(j + n, m - k) * rrho**(j + n + 1))
						cL4(icoeff, ib4i4) = cL4(icoeff, ib4i4) + rtmp * ctmp

							enddo

						enddo

						if(k.gt.0) then

							icoeff2 = igetIndex(j, -k)
							cL4(icoeff2, ib4i4) = conjg(cL4(icoeff, ib4i4))

						endif

					enddo

				enddo

			endif

			enddo
			enddo
			enddo

		enddo

	enddo

	return
	end subroutine
