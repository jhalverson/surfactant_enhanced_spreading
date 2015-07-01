	subroutine evalw2w2nn(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq, ib4
	integer*4 ib512, ib64to512
	integer*4 ictk, ictj, i, j, mi, ni, mj, mnj
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq

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

				fxw2(mi, ib512) = fxw2(mi, ib512) + xij * fel
				fyw2(mi, ib512) = fyw2(mi, ib512) + yij * fel
				fzw2(mi, ib512) = fzw2(mi, ib512) + zij * fel

			enddo

		enddo

	enddo

	return
	end subroutine
