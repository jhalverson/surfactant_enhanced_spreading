	subroutine evalSnnSnnn(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ib4
	integer*4 ictk, ictj, i, j, mi, ni, mj, nj, mni, mnj
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq
	real*8 flj, ft, rsi, r6, rpl

	do ib4 = 1, 64
	
		ictk = ictSnnCentral(ib4)
		do i = 1, ictk
		
			mni = SnnCentral(i, 1, ib4)
			mi  = SnnCentral(i, 3, ib4)
			ni  = SnnCentral(i, 4, ib4)
			
			rqq = rqs(ni) * xpi

			xi = rxs(mi, mni)
			yi = rys(mi, mni)
			zi = rzs(mi, mni)
			
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

				if(rij.le.rc) then

					rsi = 1.0D0 / rijsq
					r6 = rsi**3
					rpl = 48.0D0 * epss(ni, nj) * sgss6(ni, nj) * r6 * (sgss6(ni, nj) * r6 - 0.5D0)
					flj = rpl * rsi

				endif
				ft = fel + flj

				fxs(mi, mni) = fxs(mi, mni) + xij * ft
				fys(mi, mni) = fys(mi, mni) + yij * ft
				fzs(mi, mni) = fzs(mi, mni) + zij * ft

			enddo

		enddo

	enddo

	return
	end subroutine
