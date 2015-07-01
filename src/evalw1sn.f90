	subroutine evalw1sn(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ib4, ib216, ib64to216
	integer*4 ictk, ictj, i, j, m, n, mn
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq
	real*8 flj, ft, rsi, r6, rpl

	do ib4 = 1, 64

		ib216 = ib64to216(ib4)

		ictk = 3 * ictw1(ib216)
		do i = 1, ictk

			rqq = rqw(2) * xpi
			if(mod(i + 2, 3).eq.0) rqq = rqw(1) * xpi

			xi = rxw1(i, ib216)
			yi = ryw1(i, ib216)
			zi = rzw1(i, ib216)

			ictj = ictSnCentral(ib4)
			do j = 1, ictj

				mn = SnCentral(j, 1, ib4)
				m  = SnCentral(j, 3, ib4)
				n  = SnCentral(j, 4, ib4)

				xij = xi - rxs(m, mn)
				yij = yi - rys(m, mn)
				zij = zi - rzs(m, mn)

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
				
				fxs(m, mn) = fxs(m, mn) - xij * ft
				fys(m, mn) = fys(m, mn) - yij * ft
				fzs(m, mn) = fzs(m, mn) - zij * ft

			enddo

			ictj = ictSn26(ib4)
			do j = 1, ictj

				mn = Sn26(j, 1, ib4)
				m  = Sn26(j, 3, ib4)
				n  = Sn26(j, 4, ib4)

				xij = xi - rxs(m, mn)
				yij = yi - rys(m, mn)
				zij = zi - rzs(m, mn)

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
