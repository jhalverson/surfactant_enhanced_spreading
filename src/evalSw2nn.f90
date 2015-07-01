	subroutine evalSw2nn(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq, ib1000, ib64to1000
	integer*4 ib4
	integer*4 ictk, ictj, i, j, mi, ni, mj, nj, mnj
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq
	real*8 flj, ft, rsi, r6, rpl

	do ib4 = 1, 64
	
		ib1000 = ib64to1000(ib4)
		
		ictk = ictSCentral(ib4)
		do i = 1, ictk
		
			mi  = SCentral(i, 3, ib4)
			ni  = SCentral(i, 4, ib4)
			
			rqq = rqs(ni) * xpi

			xi = rxs(mi, ib1000)
			yi = rys(mi, ib1000)
			zi = rzs(mi, ib1000)
			
			ictj = ictW2nn(ib4)
			do j = 1, ictj

				mnj = W2nn(j, 1, ib4)
				mj  = W2nn(j, 3, ib4)

				xij = xi - rxw2(mj, mnj)
				yij = yi - ryw2(mj, mnj)
				zij = zi - rzw2(mj, mnj)

				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqw(2) / rij**3
				
				fxs(mi, ib1000) = fxs(mi, ib1000) + xij * fel
				fys(mi, ib1000) = fys(mi, ib1000) + yij * fel
				fzs(mi, ib1000) = fzs(mi, ib1000) + zij * fel
				
			enddo
			
		enddo

	enddo

	return
	end subroutine
