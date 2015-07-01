	subroutine buildW2n(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 i, j, k, l, m, n, ib512, ib64to512
	integer*4 ictk, jmolecule, kmid(26)

	real*8 rxtmp, rytmp, rztmp
	real*8 rctr4x, rctr4y, rctr4z

	logical bool26, boolCentral, isNearest, isCentral

	do i = 1, 64

		ictW2nCentral(i) = 0
		ictW2n26(i) = 0
		m = 0
		n = 0

		ib512 = ib64to512(i)

		kmid(1) = ib512 + 1
		kmid(2) = ib512 + 1 + 8
		kmid(3) = ib512 + 8
		kmid(4) = ib512 - 1 + 8
		kmid(5) = ib512 - 1
		kmid(6) = ib512 - 1 - 8
		kmid(7) = ib512 - 8
		kmid(8) = ib512 + 1 - 8
		do j = 1, 8
			kmid(j + 8) = kmid(j) + 64
			kmid(j + 16) = kmid(j) - 64
		enddo
		kmid(25) = ib512 + 64
		kmid(26) = ib512 - 64

		rctr4x = rcen4x(i)
		rctr4y = rcen4y(i)
		rctr4z = rcen4z(i)

		! scan 26 nearest neighbors
		do j = 1, 26

			k = kmid(j)
			ictk = 3 * ictw2(k)

			do l = 1, ictk

				rxtmp = rxw2(l, k)
				rytmp = ryw2(l, k)
				rztmp = rzw2(l, k)

				boolCentral = isCentral(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
				if(boolCentral) then

					m = m + 1
					jmolecule = (l - 1) / 3
					W2nCentral(m, 1, i) = k
					W2nCentral(m, 2, i) = jmolecule
					W2nCentral(m, 3, i) = l
					W2nCentral(m, 4, i) = l - 3 * jmolecule
					ictW2nCentral(i) = ictW2nCentral(i) + 1
					if(l - 3 * jmolecule.eq.1) write(*,*) "buildW2n type"

				endif

				bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
				if(bool26.and.(.not.boolCentral)) then

					n = n + 1
					jmolecule = (l - 1) / 3
					W2n26(n, 1, i) = k
					W2n26(n, 2, i) = jmolecule
					W2n26(n, 3, i) = l
					W2n26(n, 4, i) = l - 3 * jmolecule
					ictW2n26(i) = ictW2n26(i) + 1

				endif

			enddo

		enddo

	if(ictW2nCentral(i).gt.ictmaxW2n) write(*,*) "buildW2n ictmaxW2n"
	if(ictW2n26(i).gt.ictmaxW2n) write(*,*) "buildW2n ictmaxW2n"

	enddo

	return
	end subroutine
