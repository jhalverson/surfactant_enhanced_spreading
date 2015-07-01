	subroutine buildSforce(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ib4, ictk, i, mn, m, n, ib, ibox2
	integer*4 ib64m, ibox64mod, ib1000m, ib64to1000
	real*8 rx1000, ry1000, rz1000, rxL2, ryL2, rzL2

	jsplits = 0

	do ib4 = 1, 64

		ictk = ictSnCentral(ib4)
		do i = 1, ictk

			mn = SnCentral(i, 1, ib4)
			m  = SnCentral(i, 3, ib4)

			rx1000 = rcen1000x(mn)
			ry1000 = rcen1000y(mn)
			rz1000 = rcen1000z(mn)

			ib = ibox2(rx1000, ry1000, rz1000)

			if(ib.ne.iaq + 1) then

				rxL2 = rcen2x(ib)
				ryL2 = rcen2y(ib)
				rzL2 = rcen2z(ib)

				ib64m = ibox64mod(rx1000, ry1000, rz1000, rxL2, ryL2, rzL2)
				ib1000m = ib64to1000(ib64m)

				isforce(jsplits + 1) = ib - 1
				isforce(jsplits + 2) = ib1000m
				isforce(jsplits + 3) = m

				rsforce(jsplits + 1) = fxs(m, mn)
				rsforce(jsplits + 2) = fys(m, mn)
				rsforce(jsplits + 3) = fzs(m, mn)

				jsplits = jsplits + 3

			endif

		enddo

		ictk = ictSnnCentral(ib4)
		do i = 1, ictk

			mn = SnnCentral(i, 1, ib4)
			m  = SnnCentral(i, 3, ib4)

			rx1000 = rcen1000x(mn)
			ry1000 = rcen1000y(mn)
			rz1000 = rcen1000z(mn)

			ib = ibox2(rx1000, ry1000, rz1000)

			if(ib.ne.iaq + 1) then

				rxL2 = rcen2x(ib)
				ryL2 = rcen2y(ib)
				rzL2 = rcen2z(ib)

				ib64m = ibox64mod(rx1000, ry1000, rz1000, rxL2, ryL2, rzL2)
				ib1000m = ib64to1000(ib64m)

				isforce(jsplits + 1) = ib - 1
				isforce(jsplits + 2) = ib1000m
				isforce(jsplits + 3) = m

				rsforce(jsplits + 1) = fxs(m, mn)
				rsforce(jsplits + 2) = fys(m, mn)
				rsforce(jsplits + 3) = fzs(m, mn)

				jsplits = jsplits + 3

			endif

		enddo

	enddo

	if(jsplits.gt.ksplit) write(*,*) "buildSforce ksplit"
	!write(*,*) iaq, "S-initial", jsplits

	return
	end subroutine
