	subroutine evalw2nS(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ib4, ib1000, ib64to1000
	integer*4 ictk, ictj, i, j, mi, ni, mj, nj, mni, mnj
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq
	real*8 flj, ft, rsi, r6, rpl

	do ib4 = 1, 64

		ib1000 = ib64to1000(ib4)
		
		ictk = ictW2nCentral(ib4)
		do i = 1, ictk
		
			mni = W2nCentral(i, 1, ib4)
			mi  = W2nCentral(i, 3, ib4)
			
			rqq = rqw(2) * xpi

			xi = rxw2(mi, mni)
			yi = ryw2(mi, mni)
			zi = rzw2(mi, mni)

			ictj = ictSCentral(ib4)
			do j = 1, ictj

				mj  = SCentral(j, 3, ib4)
				nj  = SCentral(j, 4, ib4)

				xij = xi - rxs(mj, ib1000)
				yij = yi - rys(mj, ib1000)
				zij = zi - rzs(mj, ib1000)

				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqs(nj) / rij**3

				fxw2(mi, mni) = fxw2(mi, mni) + xij * fel
				fyw2(mi, mni) = fyw2(mi, mni) + yij * fel
				fzw2(mi, mni) = fzw2(mi, mni) + zij * fel

				fxs(mj, ib1000) = fxs(mj, ib1000) - xij * fel
				fys(mj, ib1000) = fys(mj, ib1000) - yij * fel
				fzs(mj, ib1000) = fzs(mj, ib1000) - zij * fel

			enddo
			
			ictj = ictS26(ib4)
			do j = 1, ictj

				mj  = S26(j, 3, ib4)
				nj  = S26(j, 4, ib4)

				xij = xi - rxs(mj, ib1000)
				yij = yi - rys(mj, ib1000)
				zij = zi - rzs(mj, ib1000)

				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqs(nj) / rij**3

				fxw2(mi, mni) = fxw2(mi, mni) + xij * fel
				fyw2(mi, mni) = fyw2(mi, mni) + yij * fel
				fzw2(mi, mni) = fzw2(mi, mni) + zij * fel

			enddo
			
		enddo
		
		ictk = ictSCentral(ib4)
		do i = 1, ictk
		
			mi  = SCentral(i, 3, ib4)
			ni  = SCentral(i, 4, ib4)
			
			rqq = rqs(ni) * xpi
			
			xi = rxs(mi, ib1000)
			yi = rys(mi, ib1000)
			zi = rzs(mi, ib1000)
			
			ictj = ictW2n26(ib4)
			do j = 1, ictj
			
				mnj = W2n26(j, 1, ib4)
				mj  = W2n26(j, 3, ib4)
				nj  = W2n26(j, 4, ib4)
				
				xij = xi - rxw2(mj, mnj)
				yij = yi - ryw2(mj, mnj)
				zij = zi - rzw2(mj, mnj)
				
				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqw(nj) / rij**3
				flj = 0.0D0

				if(rij.le.rc.and.nj.eq.1) then

					rsi = 1.0D0 / rijsq
					r6 = rsi**3
					rpl = 48.0D0 * epws(1, ni) * sgws6(1, ni) * r6 * (sgws6(1, ni) * r6 - 0.5D0)
					flj = rpl * rsi

				endif
				ft = fel + flj
				
				fxs(mi, ib1000) = fxs(mi, ib1000) + xij * ft
				fys(mi, ib1000) = fys(mi, ib1000) + yij * ft
				fzs(mi, ib1000) = fzs(mi, ib1000) + zij * ft
				
			enddo
			
		enddo

	enddo

	return
	end subroutine
