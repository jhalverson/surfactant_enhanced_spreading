	subroutine buildSn(iaq)
	
	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 i, j, k, l, m, n, ib1000, ib64to1000
	integer*4 ictk, jmolecule, kmid(26)

	real*8 rxtmp, rytmp, rztmp
	real*8 rctr4x, rctr4y, rctr4z

	logical bool26, boolCentral, isNearest, isCentral
	
	do i = 1, 64
	
		ictSnCentral(i) = 0
		ictSn26(i) = 0
		m = 0
		n = 0
		
		ib1000 = ib64to1000(i)
		
		kmid(1) = ib1000 + 1
		kmid(2) = ib1000 + 1 + 10
		kmid(3) = ib1000 + 10
		kmid(4) = ib1000 - 1 + 10
		kmid(5) = ib1000 - 1
		kmid(6) = ib1000 - 1 - 10
		kmid(7) = ib1000 - 10
		kmid(8) = ib1000 + 1 - 10
		do j = 1, 8
			kmid(j + 8) = kmid(j) + 100
			kmid(j + 16) = kmid(j) - 100
		enddo
		kmid(25) = ib1000 + 100
		kmid(26) = ib1000 - 100

		rctr4x = rcen4x(i)
		rctr4y = rcen4y(i)
		rctr4z = rcen4z(i)

		do j = 1, 26

			k = kmid(j)
			ictk = 29 * icts(k)

			do l = 1, ictk

				rxtmp = rxs(l, k)
				rytmp = rys(l, k)
				rztmp = rzs(l, k)

				boolCentral = isCentral(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
				if(boolCentral) then

					m = m + 1
					jmolecule = (l - 1) / 29
					SnCentral(m, 1, i) = k
					SnCentral(m, 2, i) = jmolecule
					SnCentral(m, 3, i) = l
					SnCentral(m, 4, i) = l - 29 * jmolecule
					ictSnCentral(i) = ictSnCentral(i) + 1

				endif

				bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
				if(bool26.and.(.not.boolCentral)) then

					n = n + 1
					jmolecule = (l - 1) / 29
					Sn26(n, 1, i) = k
					Sn26(n, 2, i) = jmolecule
					Sn26(n, 3, i) = l
					Sn26(n, 4, i) = l - 29 * jmolecule
					ictSn26(i) = ictSn26(i) + 1

				endif

			enddo

		enddo

		if(ictSnCentral(i).gt.ictmaxSn) write(*,*) "Bound exceeded in buildSn."
		if(ictSn26(i).gt.ictmaxSn) write(*,*) "Bound exceeded in buildSn."

	enddo

	return
	end subroutine
