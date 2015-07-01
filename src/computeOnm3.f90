	subroutine computeOnm3(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 jmid(8), jmidk(8)
	integer*4 i, j, k, n, m, jj, kk
	integer*4 icoeff, icoeff2, iIkm1, igetIndex, jmid1
	integer*4 ib216, ib8to216, ib512, ib64to512
	real*8 rctr3x, rctr3y, rctr3z, rrhopc, ralphapc, rbetapc
	real*8 rxpc, rypc, rzpc
	real*8 rAnmjk, rAnm, rtmp
	complex*16 ctmp, Ylm

	cO3 = cmplx(0.0D0, 0.0D0)

	jmid(1) = 1
	jmid(2) = 3
	jmid(3) = 9
	jmid(4) = 11
	jmid(5) = 33
	jmid(6) = 35
	jmid(7) = 41
	jmid(8) = 43

	! loop over the 8 L3 boxes
	do kk = 1, 8

		ib216 = ib8to216(kk)
		jmid1 = jmid(kk)

		jmidk(1) = jmid1
		jmidk(2) = jmid1 + 1
		jmidk(3) = jmid1 + 4
		jmidk(4) = jmid1 + 5
		jmidk(5) = jmid1 + 16
		jmidk(6) = jmid1 + 17
		jmidk(7) = jmid1 + 20
		jmidk(8) = jmid1 + 21

		rctr3x = rcen3x(kk)
		rctr3y = rcen3y(kk)
		rctr3z = rcen3z(kk)

		! loop over the 8 L4 boxes of L3 box kk
		do jj = 1, 8

			i = jmidk(jj)
			ib512 = ib64to512(i)
			rxpc = -rctr3x + rcen4x(i)
			rypc = -rctr3y + rcen4y(i)
			rzpc = -rctr3z + rcen4z(i)

			rrhopc = sqrt(rxpc**2 + rypc**2 + rzpc**2)
			ralphapc = acos(rzpc / rrhopc)
			rbetapc = atan2(rypc, rxpc)

			!if(iaq.eq.13.and.jj.eq.1.and.kk.eq.1) write(*,*) rrhopc, ralphapc, rbetapc

			do j = 0, p

				do k = 0, j

					icoeff = igetIndex(j, k)
					rAnmjk = rAnm(j, k)

					do n = 0, j

						do m = -n, n

							if(j - n - abs(k - m).ge.0) then

								ctmp = cO4(igetIndex(j - n, k - m), ib512) * Ylm(n, -m, ralphapc, rbetapc)
								rtmp = iIkm1(k, m) * rAnm(n, m) * rAnm(j - n, k - m) * rrhopc**n / rAnmjk
	!if(iaq.eq.13.and.kk.eq.1) write(*,'(i5,i5,f10.3,f10.3,f10.3,f10.3)') n, -m, cO4(igetIndex(j-n,k-m),ib512), ralphapc, rbetapc
								cO3(icoeff, ib216) = cO3(icoeff, ib216) + rtmp * ctmp

							endif

						enddo

					enddo

					if(k.gt.0) then

						icoeff2 = igetIndex(j, -k)
						cO3(icoeff2, ib216) = conjg(cO3(icoeff, ib216))

					endif

				enddo

			enddo

		enddo

	enddo

	return
	end subroutine
