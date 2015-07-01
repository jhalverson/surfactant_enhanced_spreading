	subroutine evalw2w2n(iaq)

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

			ictj = ictW2nCentral(ib4)
			do j = 1, ictj

				mnj = W2nCentral(j, 1, ib4)
				mj  = W2nCentral(j, 3, ib4)

				xij = xi - rxw2(mj, mnj)
				yij = yi - ryw2(mj, mnj)
				zij = zi - rzw2(mj, mnj)

				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqw(2) / rij**3

				fxw2(mi, ib512) = fxw2(mi, ib512) + xij * fel
				fyw2(mi, ib512) = fyw2(mi, ib512) + yij * fel
				fzw2(mi, ib512) = fzw2(mi, ib512) + zij * fel

				fxw2(mj, mnj) = fxw2(mj, mnj) - xij * fel
				fyw2(mj, mnj) = fyw2(mj, mnj) - yij * fel
				fzw2(mj, mnj) = fzw2(mj, mnj) - zij * fel

			enddo

			ictj = ictW2n26(ib4)
			do j = 1, ictj

				mnj = W2n26(j, 1, ib4)
				mj  = W2n26(j, 3, ib4)
				nj  = W2n26(j, 4, ib4)

				xij = xi - rxw2(mj, mnj)
				yij = yi - ryw2(mj, mnj)
				zij = zi - rzw2(mj, mnj)

				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqw(nj) / rij**3
				flj = 0.0D0

				if(rij.le.rc.and.ni.eq.1.and.nj.eq.1) then

					rsi = 1.0D0 / rijsq
					r6 = rsi**3
					rpl = 48.0D0 * r6 * (r6 - 0.5D0)
					flj = rpl * rsi

				endif
				ft = fel + flj

				fxw2(mi, ib512) = fxw2(mi, ib512) + xij * ft
				fyw2(mi, ib512) = fyw2(mi, ib512) + yij * ft
				fzw2(mi, ib512) = fzw2(mi, ib512) + zij * ft

			enddo

		enddo
		
		ictk = ictW2nCentral(ib4)
		do i = 1, ictk
		
			mni = W2nCentral(i, 1, ib4)
			mi  = W2nCentral(i, 3, ib4)
			
			rqq = rqw(2) * xpi

			xi = rxw2(mi, mni)
			yi = ryw2(mi, mni)
			zi = rzw2(mi, mni)

			ictj = ictW226(ib4)
			do j = 1, ictj

				mj  = W226(j, 3, ib4)
				nj  = W226(j, 4, ib4)

				xij = xi - rxw2(mj, ib512)
				yij = yi - ryw2(mj, ib512)
				zij = zi - rzw2(mj, ib512)

				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqw(nj) / rij**3

				fxw2(mi, mni) = fxw2(mi, mni) + xij * fel
				fyw2(mi, mni) = fyw2(mi, mni) + yij * fel
				fzw2(mi, mni) = fzw2(mi, mni) + zij * fel

			enddo
			
		enddo

	enddo

	return
	end subroutine
