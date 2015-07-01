	subroutine init_accel(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq, j, n, ictk
	integer*4 ib4, ib216, ib512, ib1000
	integer*4 ib64to216, ib64to512, ib64to1000

	do ib4 = 1, 64

		ib216 = ib64to216(ib4)
		ib512 = ib64to512(ib4)
		ib1000 = ib64to1000(ib4)

		ictk = 3 * ictw1(ib216)
		do j = 1, ictk

			n = j - 3 * ((j - 1) / 3)
			axw1(j, ib4) = fxw1(j, ib216) * dtsq2 / emw(n)
			ayw1(j, ib4) = fyw1(j, ib216) * dtsq2 / emw(n)
			azw1(j, ib4) = fzw1(j, ib216) * dtsq2 / emw(n)

		enddo

		ictk = 3 * ictw2(ib512)
		do j = 1, ictk

			n = j - 3 * ((j - 1) / 3)
			axw2(j, ib4) = fxw2(j, ib512) * dtsq2 / emw(n)
			ayw2(j, ib4) = fyw2(j, ib512) * dtsq2 / emw(n)
			azw2(j, ib4) = fzw2(j, ib512) * dtsq2 / emw(n)

		enddo

		ictk = 29 * icts(ib1000)
		do j = 1, ictk

			n = j - 29 * ((j - 1) / 29)
			axs(j, ib4) = fxs(j, ib1000) * dtsq2 / ems(n)
			ays(j, ib4) = fys(j, ib1000) * dtsq2 / ems(n)
			azs(j, ib4) = fzs(j, ib1000) * dtsq2 / ems(n)

		enddo

	enddo

	return
	end subroutine
