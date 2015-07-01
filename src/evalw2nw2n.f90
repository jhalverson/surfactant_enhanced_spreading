	subroutine evalw2nw2n(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq, ib4
	integer*4 ictk, ictj, i, j, mi, mj, nj, mni, mnj, mli, mlj
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq

	do ib4 = 1, 64
	
		ictk = ictW2nCentral(ib4)
		do i = 1, ictk - 1

			mni = W2nCentral(i, 1, ib4)
			mli = W2nCentral(i, 2, ib4)
			mi  = W2nCentral(i, 3, ib4)

			rqq = rqw(2) * xpi

			xi = rxw2(mi, mni)
			yi = ryw2(mi, mni)
			zi = rzw2(mi, mni)

			do j = i + 1, ictk

				mnj = W2nCentral(j, 1, ib4)
				mlj = W2nCentral(j, 2, ib4)
				mj  = W2nCentral(j, 3, ib4)

				if(mni.ne.mnj.or.mli.ne.mlj) then

					xij = xi - rxw2(mj, mnj)
					yij = yi - ryw2(mj, mnj)
					zij = zi - rzw2(mj, mnj)

					rijsq = xij**2 + yij**2 + zij**2
					rij = sqrt(rijsq)

					fel = rqq * rqw(2) / rij**3

					fxw2(mi, mni) = fxw2(mi, mni) + xij * fel
					fyw2(mi, mni) = fyw2(mi, mni) + yij * fel
					fzw2(mi, mni) = fzw2(mi, mni) + zij * fel

					fxw2(mj, mnj) = fxw2(mj, mnj) - xij * fel
					fyw2(mj, mnj) = fyw2(mj, mnj) - yij * fel
					fzw2(mj, mnj) = fzw2(mj, mnj) - zij * fel
				
				endif

			enddo
			
		enddo
		
		do i = 1, ictk

			mni = W2nCentral(i, 1, ib4)
			mli = W2nCentral(i, 2, ib4)
			mi  = W2nCentral(i, 3, ib4)

			rqq = rqw(2) * xpi

			xi = rxw2(mi, mni)
			yi = ryw2(mi, mni)
			zi = rzw2(mi, mni)

			ictj = ictW2n26(ib4)
			do j = 1, ictj

				mnj = W2n26(j, 1, ib4)
				mlj = W2n26(j, 2, ib4)
				mj  = W2n26(j, 3, ib4)
				nj  = W2n26(j, 4, ib4)

				if(mni.ne.mnj.or.mli.ne.mlj) then

					xij = xi - rxw2(mj, mnj)
					yij = yi - ryw2(mj, mnj)
					zij = zi - rzw2(mj, mnj)

					rijsq = xij**2 + yij**2 + zij**2
					rij = sqrt(rijsq)

					fel = rqq * rqw(nj) / rij**3

					fxw2(mi, mni) = fxw2(mi, mni) + xij * fel
					fyw2(mi, mni) = fyw2(mi, mni) + yij * fel
					fzw2(mi, mni) = fzw2(mi, mni) + zij * fel

				endif

			enddo

		enddo

	enddo

	return
	end subroutine
