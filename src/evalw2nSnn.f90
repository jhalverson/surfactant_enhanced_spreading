	subroutine evalw2nSnn(iaq)

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

			ictj = ictSnnCentral(ib4)
			do j = 1, ictj

				mnj = SnnCentral(j, 1, ib4)
				mj  = SnnCentral(j, 3, ib4)
				nj  = SnnCentral(j, 4, ib4)

				xij = xi - rxs(mj, mnj)
				yij = yi - rys(mj, mnj)
				zij = zi - rzs(mj, mnj)

				rijsq = xij**2 + yij**2 + zij**2
				rij = sqrt(rijsq)

				fel = rqq * rqs(nj) / rij**3

				fxw2(mi, mni) = fxw2(mi, mni) + xij * fel
				fyw2(mi, mni) = fyw2(mi, mni) + yij * fel
				fzw2(mi, mni) = fzw2(mi, mni) + zij * fel

				fxs(mj, mnj) = fxs(mj, mnj) - xij * fel
				fys(mj, mnj) = fys(mj, mnj) - yij * fel
				fzs(mj, mnj) = fzs(mj, mnj) - zij * fel

			enddo

			ictj = ictSnn26(ib4)
			do j = 1, ictj

				mnj = Snn26(j, 1, ib4)
				mj  = Snn26(j, 3, ib4)
				nj  = Snn26(j, 4, ib4)

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

		ictk = ictSnnCentral(ib4)
		do i = 1, ictk
		
			mni = SnnCentral(i, 1, ib4)
			mi  = SnnCentral(i, 3, ib4)
			ni  = SnnCentral(i, 4, ib4)
			
			rqq = rqs(ni) * xpi
			
			xi = rxs(mi, mni)
			yi = rys(mi, mni)
			zi = rzs(mi, mni)
			
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

				if(rij.le.rc.and.nj.eq.1) then

					rsi = 1.0D0 / rijsq
					r6 = rsi**3
					rpl = 48.0D0 * epws(1, ni) * sgws6(1, ni) * r6 * (sgws6(1, ni) * r6 - 0.5D0)
					flj = rpl * rsi

				endif
				ft = fel + flj
				
				fxs(mi, mni) = fxs(mi, mni) + xij * ft
				fys(mi, mni) = fys(mi, mni) + yij * ft
				fzs(mi, mni) = fzs(mi, mni) + zij * ft
				
			enddo
			
		enddo
	
	enddo

	return
	end subroutine
