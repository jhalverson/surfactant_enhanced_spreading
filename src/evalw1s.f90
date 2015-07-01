	subroutine evalw1s(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ib4, ib216, ib64to216, ib1000, ib64to1000
	integer*4 ictk, ictj, i, j, m, n
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq
	real*8 flj, ft, rsi, r6, rpl

	do ib4 = 1, 64

		ib216 = ib64to216(ib4)
		ib1000 = ib64to1000(ib4)

		ictk = 3 * ictw1(ib216)
		do i = 1, ictk

			rqq = rqw(2) * xpi
			if(mod(i + 2, 3).eq.0) rqq = rqw(1) * xpi

			xi = rxw1(i, ib216)
			yi = ryw1(i, ib216)
			zi = rzw1(i, ib216)

			ictj = ictSCentral(ib4)
			do j = 1, ictj

				m = SCentral(j, 3, ib4)
				n = SCentral(j, 4, ib4)

				xij = xi - rxs(m, ib1000)
				yij = yi - rys(m, ib1000)
				zij = zi - rzs(m, ib1000)

				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqs(n) / rij**3
				flj = 0.0D0

				if(rij.le.rc.and.mod(i + 2, 3).eq.0) then

					rsi = 1.0D0 / rijsq
					r6 = rsi**3
					rpl = 48.0D0 * epws(1, n) * sgws6(1, n) * r6 * (sgws6(1, n) * r6 - 0.5D0)
					flj = rpl * rsi

				endif
				ft = fel + flj

				fxw1(i, ib216) = fxw1(i, ib216) + xij * ft
				fyw1(i, ib216) = fyw1(i, ib216) + yij * ft
				fzw1(i, ib216) = fzw1(i, ib216) + zij * ft
				
				fxs(m, ib1000) = fxs(m, ib1000) - xij * ft
				fys(m, ib1000) = fys(m, ib1000) - yij * ft
				fzs(m, ib1000) = fzs(m, ib1000) - zij * ft

			enddo

			ictj = ictS26(ib4)
			do j = 1, ictj

				m = S26(j, 3, ib4)
				n = S26(j, 4, ib4)

				xij = xi - rxs(m, ib1000)
				yij = yi - rys(m, ib1000)
				zij = zi - rzs(m, ib1000)

				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqs(n) / rij**3
				flj = 0.0D0

				if(rij.le.rc.and.mod(i + 2, 3).eq.0) then

					rsi = 1.0D0 / rijsq
					r6 = rsi**3
					rpl = 48.0D0 * epws(1, n) * sgws6(1, n) * r6 * (sgws6(1, n) * r6 - 0.5D0)
					flj = rpl * rsi

				endif
				ft = fel + flj

				fxw1(i, ib216) = fxw1(i, ib216) + xij * ft
				fyw1(i, ib216) = fyw1(i, ib216) + yij * ft
				fzw1(i, ib216) = fzw1(i, ib216) + zij * ft

			enddo

		enddo

	enddo

	return
	end subroutine
