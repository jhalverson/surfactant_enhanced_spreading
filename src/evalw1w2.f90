	subroutine evalw1w2(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ib4, ib216, ib64to216, ib512, ib64to512
	integer*4 ictk, ictj, i, j, m, n
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq
	real*8 flj, ft, rsi, r6, rpl

	do ib4 = 1, 64

		ib216 = ib64to216(ib4)
		ib512 = ib64to512(ib4)

		ictk = 3 * ictw1(ib216)
		do i = 1, ictk

			rqq = rqw(2) * xpi
			if(mod(i + 2, 3).eq.0) rqq = rqw(1) * xpi

			xi = rxw1(i, ib216)
			yi = ryw1(i, ib216)
			zi = rzw1(i, ib216)

			ictj = ictW2Central(ib4)
			do j = 1, ictj

				m = W2Central(j, 3, ib4)
				n = W2Central(j, 4, ib4)

				xij = xi - rxw2(m, ib512)
				yij = yi - ryw2(m, ib512)
				zij = zi - rzw2(m, ib512)

				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqw(n) / rij**3
				flj = 0.0D0

				if(rij.le.rc.and.mod(i + 2, 3).eq.0.and.n.eq.1) then

					rsi = 1.0D0 / rijsq
					r6 = rsi**3
					rpl = 48.0D0 * r6 * (r6 - 0.5D0)
					flj = rpl * rsi

				endif
				ft = fel + flj

				fxw1(i, ib216) = fxw1(i, ib216) + xij * ft
				fyw1(i, ib216) = fyw1(i, ib216) + yij * ft
				fzw1(i, ib216) = fzw1(i, ib216) + zij * ft

				fxw2(m, ib512) = fxw2(m, ib512) - xij * ft
				fyw2(m, ib512) = fyw2(m, ib512) - yij * ft
				fzw2(m, ib512) = fzw2(m, ib512) - zij * ft

			enddo

			ictj = ictW226(ib4)
			do j = 1, ictj

				m = W226(j, 3, ib4)

				xij = xi - rxw2(m, ib512)
				yij = yi - ryw2(m, ib512)
				zij = zi - rzw2(m, ib512)

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
