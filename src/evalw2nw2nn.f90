	subroutine evalw2nw2nn(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq, ib4
	integer*4 ictk, ictj, i, j, mi, mj, mni, mnj
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq

	do ib4 = 1, 64
	
		ictk = ictW2nCentral(ib4)
		do i = 1, ictk
		
			mni = W2nCentral(i, 1, ib4)
			mi  = W2nCentral(i, 3, ib4)

			rqq = rqw(2) * xpi

			xi = rxw2(mi, mni)
			yi = ryw2(mi, mni)
			zi = rzw2(mi, mni)

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

				fxw2(mi, mni) = fxw2(mi, mni) + xij * fel
				fyw2(mi, mni) = fyw2(mi, mni) + yij * fel
				fzw2(mi, mni) = fzw2(mi, mni) + zij * fel

			enddo
			
		enddo

	enddo

	return
	end subroutine
