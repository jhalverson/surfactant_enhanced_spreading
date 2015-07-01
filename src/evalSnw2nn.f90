	subroutine evalSnw2nn(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ib4
	integer*4 ictk, ictj, i, j, mi, ni, mj, nj, mni, mnj
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq
	real*8 flj, ft, rsi, r6, rpl

	do ib4 = 1, 64
	
		ictk = ictSnCentral(ib4)
		do i = 1, ictk
		
			mni = SnCentral(i, 1, ib4)
			mi  = SnCentral(i, 3, ib4)
			ni  = SnCentral(i, 4, ib4)
			
			rqq = rqs(ni) * xpi

			xi = rxs(mi, mni)
			yi = rys(mi, mni)
			zi = rzs(mi, mni)
			
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

				fxs(mi, mni) = fxs(mi, mni) + xij * fel
				fys(mi, mni) = fys(mi, mni) + yij * fel
				fzs(mi, mni) = fzs(mi, mni) + zij * fel

			enddo
			
		enddo

	enddo

	return
	end subroutine
