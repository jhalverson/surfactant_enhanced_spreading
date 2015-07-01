	subroutine buildS(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 i, j, k, l, m, n, ib1000, ib64to1000
	integer*4 ictk, jmolecule

	real*8 rxtmp, rytmp, rztmp
	real*8 rctr4x, rctr4y, rctr4z

	logical bool26, boolCentral, isnearest, isCentral

	do i = 1, 64

		ictSCentral(i) = 0
		ictS26(i) = 0
		m = 0
		n = 0

		ib1000 = ib64to1000(i)

		rctr4x = rcen4x(i)
		rctr4y = rcen4y(i)
		rctr4z = rcen4z(i)

		k = ib1000
		ictk = 29 * icts(k)

		do l = 1, ictk

			rxtmp = rxs(l, k)
			rytmp = rys(l, k)
			rztmp = rzs(l, k)

			boolCentral = isCentral(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
			if(boolCentral) then

				m = m + 1
				jmolecule = (l - 1) / 29
				SCentral(m, 1, i) = k
				SCentral(m, 2, i) = jmolecule
				SCentral(m, 3, i) = l
				SCentral(m, 4, i) = l - 29 * jmolecule
				ictSCentral(i) = ictSCentral(i) + 1

			endif

			bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
			if(bool26.and.(.not.boolCentral)) then

				n = n + 1
				jmolecule = (l - 1) / 29
				S26(n, 1, i) = k
				S26(n, 2, i) = jmolecule
				S26(n, 3, i) = l
				S26(n, 4, i) = l - 29 * jmolecule
				ictS26(i) = ictS26(i) + 1
				!write(*,'(f8.4,f8.4,f8.4,f8.4,f8.4,f8.4)') rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z

			endif

		enddo
		!write(*,*) ictk, m, n
	if(ictSCentral(i).gt.ictmaxS) write(*,*) "Bound exceeded in buildS."
	if(ictS26(i).gt.ictmaxS) write(*,*) "Bound exceeded in buildS."
	!if(iaq.eq.1) write(*,*) i, ib1000, ictSCentral(i), ictS26(i)

	enddo

	return
	end subroutine
