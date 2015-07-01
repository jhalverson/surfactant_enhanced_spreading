	subroutine computeLjk2(iaq, ibq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"

	integer*4 iaq, ibq(3)
	integer*4 ii, jj, kk, j, k, l, m, n
	integer*4 icoeff, icoeff2, iIkm2, igetIndex
	integer*4 ibq3, ibq2, ibq1
	real*8 rctr2x, rctr2y, rctr2z
	real*8 rcendiffx, rcendiffy, rcendiffz
	real*8 rrho, ralpha, rbeta
	real*8 rAnmjk, rAnm, rtmp
	complex*16 ctmp, Ylm

	cL2 = cmplx(0.0D0, 0.0D0)

	ibq3 = ibq(3) + 1
	ibq2 = ibq(2) + 1
	ibq1 = ibq(1) + 1

	do kk = 1, ibdz2
	do jj = 1, ibdy2
	do ii = 1, ibdx2

	if(abs(ii - ibq3).gt.1.or.abs(jj - ibq2).gt.1.or.abs(kk - ibq1).gt.1) then

		l = ii + (jj - 1) * ibdx2 + (kk - 1) * ibdx2 * ibdy2
		rcendiffx = rcen2x(l) - rcen2xtmp
		rcendiffy = rcen2y(l) - rcen2ytmp
		rcendiffz = rcen2z(l) - rcen2ztmp

		rrho = sqrt(rcendiffx**2 + rcendiffy**2 + rcendiffz**2)
		ralpha = acos(rcendiffz / rrho)
		rbeta = atan2(rcendiffy, rcendiffx)

		do j = 0, p

			do k = 0, j

				icoeff = igetIndex(j, k)
				rAnmjk = rAnm(j, k)

				do n = 0, p

					do m = -n, n

						ctmp = cO2all(igetIndex(n, m) + lp * (l - 1)) * Ylm(j + n, m - k, ralpha, rbeta)
						rtmp = iIkm2(k, m) * rAnmjk * rAnm(n, m) / ((-1)**n * rAnm(j + n, m - k) * rrho**(j + n + 1))
						cL2(icoeff) = cL2(icoeff) + rtmp * ctmp

					enddo

				enddo

				if(k.gt.0) then

					icoeff2 = igetIndex(j, -k)
					cL2(icoeff2) = conjg(cL2(icoeff))

				endif

			enddo

		enddo

	endif

	enddo
	enddo
	enddo

	return
	end subroutine
