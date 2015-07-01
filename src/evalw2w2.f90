	subroutine evalw2w2(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ib4, ib512, ib64to512
	integer*4 ictk, ictj, i, j
	integer*4 mli, mlj, mi, mj, ni, nj
	real*8 xi, yi, zi, xij, yij, zij
	real*8 rijsq, rij, fel, rqq
	real*8 flj, ft, rsi, r6, rpl

	do ib4 = 1, 64

		ib512 = ib64to512(ib4)

		ictk = ictW2Central(ib4)
		do i = 1, ictk - 1

			mli = W2Central(i, 2, ib4)
			mi  = W2Central(i, 3, ib4)
			ni  = W2Central(i, 4, ib4)

			rqq = rqw(ni) * xpi

			xi = rxw2(mi, ib512)
			yi = ryw2(mi, ib512)
			zi = rzw2(mi, ib512)

			do j = i + 1, ictk

				mlj = W2Central(j, 2, ib4)
				mj  = W2Central(j, 3, ib4)
				nj  = W2Central(j, 4, ib4)

				if(mli.ne.mlj) then

					xij = xi - rxw2(mj, ib512)
					yij = yi - ryw2(mj, ib512)
					zij = zi - rzw2(mj, ib512)

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

					fxw2(mj, ib512) = fxw2(mj, ib512) - xij * ft
					fyw2(mj, ib512) = fyw2(mj, ib512) - yij * ft
					fzw2(mj, ib512) = fzw2(mj, ib512) - zij * ft

				endif

			enddo

		enddo

		do i = 1, ictk

			mli = W2Central(i, 2, ib4)
			mi  = W2Central(i, 3, ib4)
			ni  = W2Central(i, 4, ib4)

			rqq = rqw(ni) * xpi

			xi = rxw2(mi, ib512)
			yi = ryw2(mi, ib512)
			zi = rzw2(mi, ib512)

			ictj = ictW226(ib4)
			do j = 1, ictj

				mlj = W226(j, 2, ib4)
				mj  = W226(j, 3, ib4)

				if(mli.ne.mlj) then

					xij = xi - rxw2(mj, ib512)
					yij = yi - ryw2(mj, ib512)
					zij = zi - rzw2(mj, ib512)

					rijsq = xij**2 + yij**2 + zij**2
					rij = sqrt(rijsq)

					fel = rqq * rqw(2) / rij**3

					fxw2(mi, ib512) = fxw2(mi, ib512) + xij * fel
					fyw2(mi, ib512) = fyw2(mi, ib512) + yij * fel
					fzw2(mi, ib512) = fzw2(mi, ib512) + zij * fel

				endif

			enddo

		enddo

	enddo

	return
	end subroutine
