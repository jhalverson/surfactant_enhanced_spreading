	subroutine buildW2(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 i, j, k, l, m, n, ib512, ib64to512
	integer*4 ictk, jmolecule
	integer*4 ibox64

	real*8 rxtmp, rytmp, rztmp
	real*8 rctr4x, rctr4y, rctr4z

	logical bool26, boolCentral, isNearest, isCentral

	do i = 1, 64

		ictW2Central(i) = 0
		ictW226(i) = 0
		m = 0
		n = 0

		ib512 = ib64to512(i)

		rctr4x = rcen4x(i)
		rctr4y = rcen4y(i)
		rctr4z = rcen4z(i)

		k = ib512
		ictk = 3 * ictw2(k)

		do l = 1, ictk

			rxtmp = rxw2(l, k)
			rytmp = ryw2(l, k)
			rztmp = rzw2(l, k)

			boolCentral = isCentral(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
			if(boolCentral) then

				m = m + 1
				jmolecule = (l - 1) / 3
				W2Central(m, 1, i) = k
				W2Central(m, 2, i) = jmolecule
				W2Central(m, 3, i) = l
				W2Central(m, 4, i) = l - 3 * jmolecule
				ictW2Central(i) = ictW2Central(i) + 1

			endif

			bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
			if(bool26.and.(.not.boolCentral)) then

				n = n + 1 
				jmolecule = (l - 1) / 3
				W226(n, 1, i) = k
				W226(n, 2, i) = jmolecule
				W226(n, 3, i) = l
				W226(n, 4, i) = l - 3 * jmolecule
				ictW226(i) = ictW226(i) + 1
				if(l - 3 * jmolecule.eq.1) write(*,*) "buildW2 type"

			endif

		enddo

	if(ictW2Central(i).gt.ictmaxW2) write(*,*) "buildW2 ictmaxW2"
	if(ictW226(i).gt.ictmaxW2) write(*,*) "buildW2 ictmaxW2"
	if(ictW2Central(i) + ictW226(i).ne.ictk) write(*,*) "buildW2 ictk"
	
	enddo

	return
	end subroutine
