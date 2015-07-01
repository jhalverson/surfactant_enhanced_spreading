	subroutine buildW2force(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq
	integer*4 ib4, ictk, i, mn, m, n, ib, ibox2
	integer*4 ib64m, ibox64mod, ib512, ib64to512
	real*8 rx512, ry512, rz512, rxL2, ryL2, rzL2

	integer*4 jmolecule, jtype

	jsplitw2 = 0

	do ib4 = 1, 64

		ictk = ictW2nCentral(ib4)
		do i = 1, ictk

			mn = W2nCentral(i, 1, ib4)
			m  = W2nCentral(i, 3, ib4)

			rx512 = rcen512x(mn)
			ry512 = rcen512y(mn)
			rz512 = rcen512z(mn)

			ib = ibox2(rx512, ry512, rz512)

			if(ib.ne.iaq + 1) then

				rxL2 = rcen2x(ib)
				ryL2 = rcen2y(ib)
				rzL2 = rcen2z(ib)

				ib64m = ibox64mod(rx512, ry512, rz512, rxL2, ryL2, rzL2)
				ib512 = ib64to512(ib64m)

				iw2force(jsplitw2 + 1) = ib - 1
				iw2force(jsplitw2 + 2) = ib512
				iw2force(jsplitw2 + 3) = m
				jmolecule = (m - 1) / 3
				jtype = m - 3 * jmolecule
				!write(*,*) jtype, jmolecule, m, mn
				if(jtype.ne.2.and.jtype.ne.3) write(*,*) "buildW2force jtype"

				rw2force(jsplitw2 + 1) = fxw2(m, mn)
				rw2force(jsplitw2 + 2) = fyw2(m, mn)
				rw2force(jsplitw2 + 3) = fzw2(m, mn)

				jsplitw2 = jsplitw2 + 3

			endif

		enddo

	enddo

	if(jsplitw2.gt.ksplit) write(*,*) "buildW2force ksplit"
	!write(*,*) iaq, "w2-initial", jsplitw2

	return
	end subroutine
