	subroutine evalw2sn(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ib4, ib512, ib64to512
	integer*4 ictk, ictj, i, j, mi, ni, mni, mj, nj, mnj
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq
	real*8 flj, ft, rsi, r6, rpl

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

			ictj = ictSnCentral(ib4)
			do j = 1, ictj

				mnj = SnCentral(j, 1, ib4)
				mj  = SnCentral(j, 3, ib4)
				nj  = SnCentral(j, 4, ib4)

				xij = xi - rxs(mj, mnj)
				yij = yi - rys(mj, mnj)
				zij = zi - rzs(mj, mnj)

				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqs(nj) / rij**3
				flj = 0.0D0

				if(rij.le.rc.and.ni.eq.1) then

					rsi = 1.0D0 / rijsq
					r6 = rsi**3
					rpl = 48.0D0 * epws(1, nj) * sgws6(1, nj) * r6 * (sgws6(1, nj) * r6 - 0.5D0)
					flj = rpl * rsi

				endif
				ft = fel + flj

				fxw2(mi, ib512) = fxw2(mi, ib512) + xij * ft
				fyw2(mi, ib512) = fyw2(mi, ib512) + yij * ft
				fzw2(mi, ib512) = fzw2(mi, ib512) + zij * ft
				
				fxs(mj, mnj) = fxs(mj, mnj) - xij * ft
				fys(mj, mnj) = fys(mj, mnj) - yij * ft
				fzs(mj, mnj) = fzs(mj, mnj) - zij * ft

			enddo

			ictj = ictSn26(ib4)
			do j = 1, ictj

				mnj = Sn26(j, 1, ib4)
				mj  = Sn26(j, 3, ib4)
				nj  = Sn26(j, 4, ib4)

				xij = xi - rxs(mj, mnj)
				yij = yi - rys(mj, mnj)
				zij = zi - rzs(mj, mnj)

				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqs(nj) / rij**3
				flj = 0.0D0

				if(rij.le.rc.and.ni.eq.1) then

					rsi = 1.0D0 / rijsq
					r6 = rsi**3
					rpl = 48.0D0 * epws(1, nj) * sgws6(1, nj) * r6 * (sgws6(1, nj) * r6 - 0.5D0)
					flj = rpl * rsi

				endif
				ft = fel + flj

				fxw2(mi, ib512) = fxw2(mi, ib512) + xij * ft
				fyw2(mi, ib512) = fyw2(mi, ib512) + yij * ft
				fzw2(mi, ib512) = fzw2(mi, ib512) + zij * ft

			enddo

		enddo

		ictk = ictSnCentral(ib4)
		do i = 1, ictk
		
			mni = SnCentral(i, 1, ib4)
			mi  = SnCentral(i, 3, ib4)
			ni  = SnCentral(i, 4, ib4)
			
			rqq = rqs(ni) * xpi
			
			xi = rxs(mi, mni)
			yi = rys(mi, mni)
			zi = rzs(mi, mni)
			
			ictj = ictW226(ib4)
			do j = 1, ictj
			
				mj = W226(j, 3, ib4)
				
				xij = xi - rxw2(mj, ib512)
				yij = yi - ryw2(mj, ib512)
				zij = zi - rzw2(mj, ib512)
				
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
