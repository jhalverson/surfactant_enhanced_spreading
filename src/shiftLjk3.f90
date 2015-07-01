	subroutine shiftLjk3(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 kk, j, k, n, m
	integer*4 icoeff, icoeff2, iIkm3, igetIndex
	integer*4 key4(8), ib4(8), i, i3, i4, key41
	real*8 rxpc, rypc, rzpc
	real*8 rrhopc, ralphapc, rbetapc
	real*8 rctr3x, rctr3y, rctr3z
	real*8 rAnmjk, rAnm, rtmp
	complex*16 ctmp, Ylm

	cL4 = cmplx(0.0D0, 0.0D0)

	key4(1) = 1
	key4(2) = 3
	key4(3) = 9
	key4(4) = 11
	key4(5) = 33
	key4(6) = 35
	key4(7) = 41
	key4(8) = 43

	do i3 = 1, 8

		key41 = key4(i3)

		ib4(1) = key41
		ib4(2) = key41 + 1
		ib4(3) = key41 + 4
		ib4(4) = key41 + 5
		ib4(5) = key41 + 16
		ib4(6) = key41 + 17
		ib4(7) = key41 + 20
		ib4(8) = key41 + 21

		rctr3x = rcen3x(i3)
		rctr3y = rcen3y(i3)
		rctr3z = rcen3z(i3)

		do i4 = 1, 8

			i = ib4(i4)
			rxpc = rctr3x - rcen4x(i)
			rypc = rctr3y - rcen4y(i)
			rzpc = rctr3z - rcen4z(i)

			rrhopc = sqrt(rxpc**2 + rypc**2 + rzpc**2)
			ralphapc = acos(rzpc / rrhopc)
			rbetapc = atan2(rypc, rxpc)

			do j = 0, p

				do k = 0, j

					icoeff = igetIndex(j, k)
					rAnmjk = rAnm(j, k)

					do n = 0, p

						do m = -n, n

						if(n - j - abs(m - k).ge.0) then

						ctmp = cL3(igetIndex(n, m), i3) * Ylm(n - j, m - k, ralphapc, rbetapc)
						rtmp = iIkm3(k, m) * rAnm(n - j, m - k) * rAnmjk * rrhopc**(n - j) / ((-1)**(n + j) * rAnm(n, m))
						cL4(icoeff, i) = cL4(icoeff, i) + rtmp * ctmp

						endif

						enddo

					enddo

					if(k.gt.0) then

						icoeff2 = igetIndex(j, -k)
						cL4(icoeff2, i) = conjg(cL4(icoeff, i))

					endif

				enddo

			enddo

		enddo

	enddo

	return
	end subroutine
