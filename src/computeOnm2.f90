	subroutine computeOnm2(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 j, k, n, m, kk
	integer*4 icoeff, icoeff2, iIkm1, igetIndex, ib216, ib8to216
	real*8 rrhopc, ralphapc, rbetapc
	real*8 rxpc, rypc, rzpc
	real*8 rAnmjk, rAnm, rtmp
	complex*16 ctmp, Ylm

	cO2 = cmplx(0.0D0, 0.0D0)

	do kk = 1, 8

		ib216 = ib8to216(kk)

		rxpc = -rcen2xtmp + rcen3x(kk)
		rypc = -rcen2ytmp + rcen3y(kk)
		rzpc = -rcen2ztmp + rcen3z(kk)

		rrhopc = sqrt(rxpc**2 + rypc**2 + rzpc**2)
		ralphapc = acos(rzpc / rrhopc)
		rbetapc = atan2(rypc, rxpc)

		do j = 0, p

			do k = 0, j

				icoeff = igetIndex(j, k)
				rAnmjk = rAnm(j, k)

				do n = 0, j

					do m = -n, n

						if(j - n - abs(k - m).ge.0) then

							ctmp = cO3(igetIndex(j - n, k - m), ib216) * Ylm(n, -m, ralphapc, rbetapc)
							rtmp = iIkm1(k, m) * rAnm(n, m) * rAnm(j - n, k - m) * rrhopc**n / rAnmjk
							cO2(icoeff) = cO2(icoeff) + rtmp * ctmp

						endif

					enddo

				enddo

				if(k.gt.0) then

					icoeff2 = igetIndex(j, -k)
					cO2(icoeff2) = conjg(cO2(icoeff))

				endif

			enddo

		enddo

	enddo

	return
	end subroutine
