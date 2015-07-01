	subroutine evalw1w2nn(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ib4, ib216, ib64to216
	integer*4 ictk, ictj, i, j, m, mn
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq

	do ib4 = 1, 64

		ib216 = ib64to216(ib4)

		ictk = 3 * ictw1(ib216)
		do i = 1, ictk

			rqq = rqw(2) * xpi
			if(mod(i + 2, 3).eq.0) rqq = rqw(1) * xpi

			xi = rxw1(i, ib216)
			yi = ryw1(i, ib216)
			zi = rzw1(i, ib216)

			ictj = ictW2nn(ib4)
			do j = 1, ictj

				mn = W2nn(j, 1, ib4)
				m  = W2nn(j, 3, ib4)

				xij = xi - rxw2(m, mn)
				yij = yi - ryw2(m, mn)
				zij = zi - rzw2(m, mn)

				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqw(2) / rij**3

				fxw1(i, ib216) = fxw1(i, ib216) + xij * fel
				fyw1(i, ib216) = fyw1(i, ib216) + yij * fel
				fzw1(i, ib216) = fzw1(i, ib216) + zij * fel

			enddo

		enddo

	enddo

	return
	end subroutine
