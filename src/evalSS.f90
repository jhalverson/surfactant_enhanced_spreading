	subroutine evalSS(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq, ib1000, ib64to1000
	integer*4 ib4
	integer*4 ictk, ictj, i, j, mi, ni, mj, nj, mli, mlj
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq
	real*8 flj, ft, rsi, r6, rpl

	do ib4 = 1, 64
	
		ib1000 = ib64to1000(ib4)
		
		ictk = ictSCentral(ib4)
		do i = 1, ictk - 1
		
			mli = SCentral(i, 2, ib4)
			mi  = SCentral(i, 3, ib4)
			ni  = SCentral(i, 4, ib4)
			
			rqq = rqs(ni) * xpi

			xi = rxs(mi, ib1000)
			yi = rys(mi, ib1000)
			zi = rzs(mi, ib1000)
			
			do j = i + 1, ictk

				mlj = SCentral(j, 2, ib4)
				mj  = SCentral(j, 3, ib4)
				nj  = SCentral(j, 4, ib4)
				
				if(mli.ne.mlj) then

					xij = xi - rxs(mj, ib1000)
					yij = yi - rys(mj, ib1000)
					zij = zi - rzs(mj, ib1000)

					rijsq = xij**2 + yij**2 + zij**2
					rij = sqrt(rijsq)

					fel = rqq * rqs(nj) / rij**3
					flj = 0.0D0

					if(rij.le.rc) then

						rsi = 1.0D0 / rijsq
						r6 = rsi**3
						rpl = 48.0D0 * epss(ni, nj) * sgss6(ni, nj) * r6 * (sgss6(ni, nj) * r6 - 0.5D0)
						flj = rpl * rsi

					endif
					ft = fel + flj

					fxs(mi, ib1000) = fxs(mi, ib1000) + xij * ft
					fys(mi, ib1000) = fys(mi, ib1000) + yij * ft
					fzs(mi, ib1000) = fzs(mi, ib1000) + zij * ft
					
					fxs(mj, ib1000) = fxs(mj, ib1000) - xij * ft
					fys(mj, ib1000) = fys(mj, ib1000) - yij * ft
					fzs(mj, ib1000) = fzs(mj, ib1000) - zij * ft

				endif

			enddo
			
		enddo
		
		do i = 1, ictk
		
			mli = SCentral(i, 2, ib4)
			mi  = SCentral(i, 3, ib4)
			ni  = SCentral(i, 4, ib4)
			
			rqq = rqs(ni) * xpi

			xi = rxs(mi, ib1000)
			yi = rys(mi, ib1000)
			zi = rzs(mi, ib1000)
			
			ictj = ictS26(ib4)
			do j = 1, ictj

				mlj = S26(j, 2, ib4)
				mj  = S26(j, 3, ib4)
				nj  = S26(j, 4, ib4)
				
				if(mli.ne.mlj) then

					xij = xi - rxs(mj, ib1000)
					yij = yi - rys(mj, ib1000)
					zij = zi - rzs(mj, ib1000)

					rijsq = xij**2 + yij**2 + zij**2
					rij = sqrt(rijsq)

					fel = rqq * rqs(nj) / rij**3
					flj = 0.0D0

					if(rij.le.rc) then

						rsi = 1.0D0 / rijsq
						r6 = rsi**3
						rpl = 48.0D0 * epss(ni, nj) * sgss6(ni, nj) * r6 * (sgss6(ni, nj) * r6 - 0.5D0)
						flj = rpl * rsi

					endif

					ft = fel + flj

					fxs(mi, ib1000) = fxs(mi, ib1000) + xij * ft
					fys(mi, ib1000) = fys(mi, ib1000) + yij * ft
					fzs(mi, ib1000) = fzs(mi, ib1000) + zij * ft

				endif

			enddo
			
		enddo

	enddo

	return
	end subroutine
