	subroutine evalw2w1n(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ib4, ib216, ib64to216, ib512, ib64to512
	integer*4 kmid(26), ib, ib26, ictj
	integer*4 ictk, i, j, j1, j2
	integer*4 mi, ni
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq
	real*8 flj, ft, rsi, r6, rpl

	do ib4 = 1, 64

		ib216 = ib64to216(ib4)
		ib512 = ib64to512(ib4)
		
		kmid(1) = ib216 + 1
		kmid(2) = ib216 + 1 + 6
		kmid(3) = ib216 + 6
		kmid(4) = ib216 - 1 + 6
		kmid(5) = ib216 - 1
		kmid(6) = ib216 - 1 - 6
		kmid(7) = ib216 - 6
		kmid(8) = ib216 + 1 - 6
		do j = 1, 8
			kmid(j + 8) = kmid(j) + 36
			kmid(j + 16) = kmid(j) - 36
		enddo
		kmid(25) = ib216 + 36
		kmid(26) = ib216 - 36
		
		ictk = ictW2Central(ib4)
		do i = 1, ictk

			mi  = W2Central(i, 3, ib4)
			ni  = W2Central(i, 4, ib4)

			rqq = rqw(ni) * xpi

			xi = rxw2(mi, ib512)
			yi = ryw2(mi, ib512)
			zi = rzw2(mi, ib512)
			
			do ib = 1, 26
			
				ib26 = kmid(ib)

				ictj = 3 * ictw1(ib26)
				do j = 1, ictj - 2, 3

					xij = xi - rxw1(j, ib26)
					yij = yi - ryw1(j, ib26)
					zij = zi - rzw1(j, ib26)

					rijsq = xij**2 + yij**2 + zij**2
					rij = sqrt(rijsq)

					fel = rqq * rqw(1) / rij**3
					flj = 0.0D0

					if(rij.le.rc.and.ni.eq.1) then

						rsi = 1.0D0 / rijsq
						r6 = rsi**3
						rpl = 48.0D0 * r6 * (r6 - 0.5D0)
						flj = rpl * rsi

					endif
					ft = fel + flj

					fxw2(mi, ib512) = fxw2(mi, ib512) + xij * ft
					fyw2(mi, ib512) = fyw2(mi, ib512) + yij * ft
					fzw2(mi, ib512) = fzw2(mi, ib512) + zij * ft
					
					j1 = j + 1
					xij = xi - rxw1(j1, ib26)
					yij = yi - ryw1(j1, ib26)
					zij = zi - rzw1(j1, ib26)

					rijsq = xij**2 + yij**2 + zij**2
					rij = sqrt(rijsq)

					fel = rqq * rqw(2) / rij**3
					
					fxw2(mi, ib512) = fxw2(mi, ib512) + xij * fel
					fyw2(mi, ib512) = fyw2(mi, ib512) + yij * fel
					fzw2(mi, ib512) = fzw2(mi, ib512) + zij * fel
					
					j2 = j + 2
					xij = xi - rxw1(j2, ib26)
					yij = yi - ryw1(j2, ib26)
					zij = zi - rzw1(j2, ib26)

					rijsq = xij**2 + yij**2 + zij**2
					rij = sqrt(rijsq)

					fel = rqq * rqw(3) / rij**3

					fxw2(mi, ib512) = fxw2(mi, ib512) + xij * fel
					fyw2(mi, ib512) = fyw2(mi, ib512) + yij * fel
					fzw2(mi, ib512) = fzw2(mi, ib512) + zij * fel

				enddo

			enddo

		enddo

	enddo

	return
	end subroutine
