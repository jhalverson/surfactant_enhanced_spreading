	subroutine shiftLjk2(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 kk, j, k, n, m
	integer*4 icoeff, icoeff2, iIkm3, igetIndex
	real*8 rxpc, rypc, rzpc
	real*8 rrhopc, ralphapc, rbetapc
	real*8 rAnmjk, rAnm, rtmp
	complex*16 ctmp, Ylm

	cL3 = cmplx(0.0D0, 0.0D0)

	do kk = 1, 8

		rxpc = rcen2xtmp - rcen3x(kk)
		rypc = rcen2ytmp - rcen3y(kk)
		rzpc = rcen2ztmp - rcen3z(kk)

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

						ctmp = cL2(igetIndex(n, m)) * Ylm(n - j, m - k, ralphapc, rbetapc)
						rtmp = iIkm3(k, m) * rAnm(n - j, m - k) * rAnmjk * rrhopc**(n - j) / ((-1)**(n + j) * rAnm(n, m))
						cL3(icoeff, kk) = cL3(icoeff, kk) + rtmp * ctmp

						endif

					enddo

				enddo

				if(k.gt.0) then

					icoeff2 = igetIndex(j, -k)
					cL3(icoeff2, kk) = conjg(cL3(icoeff, kk))

				endif

			enddo

		enddo

	enddo

	return
	end subroutine
