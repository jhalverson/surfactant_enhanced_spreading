	subroutine evalw2nSnnn(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ib4
	integer*4 ictk, ictj, i, j, mi, ni, mj, nj, mni, mnj
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq
	real*8 flj, ft, rsi, r6, rpl

	do ib4 = 1, 64
	
		ictk = ictW2nCentral(ib4)
		do i = 1, ictk
		
			mni = W2nCentral(i, 1, ib4)
			mi  = W2nCentral(i, 3, ib4)
			
			rqq = rqw(2) * xpi

			xi = rxw2(mi, mni)
			yi = ryw2(mi, mni)
			zi = rzw2(mi, mni)

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

				fxw2(mi, mni) = fxw2(mi, mni) + xij * fel
				fyw2(mi, mni) = fyw2(mi, mni) + yij * fel
				fzw2(mi, mni) = fzw2(mi, mni) + zij * fel

			enddo

		enddo

	enddo

	return
	end subroutine
