	subroutine evalw2snnn(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ib4, ib512, ib64to512
	integer*4 ictk, ictj, i, j, mi, ni, mj, nj, mnj
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq
	real*8 flj, ft, rsi, r6, rpl

	do ib4 = 1, 64

		ib512 = ib64to512(ib4)

		ictk = ictW2Central(ib4)
		do i = 1, ictk
		
			mi  = W2Central(i, 3, ib4)
			ni  = W2Central(i, 4, ib4)
			
			rqq = rqw(ni) * xpi

			xi = rxw2(mi, ib512)
			yi = ryw2(mi, ib512)
			zi = rzw2(mi, ib512)

			ictj = ictSnnn(ib4)
			do j = 1, ictj

				mnj = Snnn(j, 1, ib4)
				mj  = Snnn(j, 3, ib4)
				nj  = Snnn(j, 4, ib4)

				xij = xi - rxs(mj, mnj)
				yij = yi - rys(mj, mnj)
				zij = zi - rzs(mj, mnj)

				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqs(nj) / rij**3
				flj = 0.0D0

				if(rij.le.rc.and.ni.eq.1) then

					rsi = 1.0D0 / rijsq
					r6 = rsi**3
					rpl = 48.0D0 * epws(1, nj) * sgws6(1, nj) * r6 * (sgws6(1, nj) * r6 - 0.5D0)
					flj = rpl * rsi

				endif
				ft = fel + flj

				fxw2(mi, ib512) = fxw2(mi, ib512) + xij * ft
				fyw2(mi, ib512) = fyw2(mi, ib512) + yij * ft
				fzw2(mi, ib512) = fzw2(mi, ib512) + zij * ft

			enddo

		enddo

	enddo

	return
	end subroutine
