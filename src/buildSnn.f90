	! check 5^3 - 3^3 = 98 boxes
	subroutine buildSnn(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 i, j, k, l, m, n, ib1000, ib64to1000
	integer*4 ibstart, ibend, ictk, jmolecule

	real*8 rxtmp, rytmp, rztmp
	real*8 rctr4x, rctr4y, rctr4z

	logical bool26, boolCentral, isnearest, isCentral

	do i = 1, 64

		ictSnnCentral(i) = 0
		ictSnn26(i) = 0
		m = 0
		n = 0

		ib1000 = ib64to1000(i)
		rctr4x = rcen4x(i)
		rctr4y = rcen4y(i)
		rctr4z = rcen4z(i)

		! lower 5 by 5
		do j = 1, 5

			ibstart = ib1000 - 200 - 22 + (j - 1) * 10
			ibend = ibstart + 4

			do k = ibstart, ibend
			!if(iaq.eq.13.and.i.eq.1) write(*,*) ib1000, k

				ictk = 29 * icts(k)

				do l = 1, ictk

					rxtmp = rxs(l, k)
					rytmp = rys(l, k)
					rztmp = rzs(l, k)

					boolCentral = isCentral(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
					if(boolCentral) then

						m = m + 1
						jmolecule = (l - 1) / 29
						SnnCentral(m, 1, i) = k
						SnnCentral(m, 2, i) = jmolecule
						SnnCentral(m, 3, i) = l
						SnnCentral(m, 4, i) = l - 29 * jmolecule
						ictSnnCentral(i) = ictSnnCentral(i) + 1

					endif

					bool26 = isnearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
					if(bool26.and.(.not.boolCentral)) then

						n = n + 1
						jmolecule = (l - 1) / 29
						Snn26(n, 1, i) = k
						Snn26(n, 2, i) = jmolecule
						Snn26(n, 3, i) = l
						Snn26(n, 4, i) = l - 29 * jmolecule
						ictSnn26(i) = ictSnn26(i) + 1

					endif

				enddo

			enddo

		enddo

		! upper 5 by 5
		do j = 1, 5

			ibstart = ib1000 + 200 - 22 + (j - 1) * 10
			ibend = ibstart + 4

			do k = ibstart, ibend
			!if(iaq.eq.13.and.i.eq.1) write(*,*) ib1000, k

				ictk = 29 * icts(k)

				do l = 1, ictk

					rxtmp = rxs(l, k)
					rytmp = rys(l, k)
					rztmp = rzs(l, k)

					boolCentral = isCentral(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
					if(boolCentral) then

						m = m + 1
						jmolecule = (l - 1) / 29
						SnnCentral(m, 1, i) = k
						SnnCentral(m, 2, i) = jmolecule
						SnnCentral(m, 3, i) = l
						SnnCentral(m, 4, i) = l - 29 * jmolecule
						ictSnnCentral(i) = ictSnnCentral(i) + 1

					endif

					bool26 = isnearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
					if(bool26.and.(.not.boolCentral)) then

						n = n + 1
						jmolecule = (l - 1) / 29
						Snn26(n, 1, i) = k
						Snn26(n, 2, i) = jmolecule
						Snn26(n, 3, i) = l
						Snn26(n, 4, i) = l - 29 * jmolecule
						ictSnn26(i) = ictSnn26(i) + 1

					endif

				enddo

			enddo

		enddo
		
		! square rings
		do j = -1, 1

			! south boxes (constant -y)
			ibstart = ib1000 + 100 * j - 20 - 2
			ibend = ibstart + 4

			do k = ibstart, ibend
			!if(iaq.eq.13.and.i.eq.1) write(*,*) ib1000, k

				ictk = 29 * icts(k)

				do l = 1, ictk
				
					rxtmp = rxs(l, k)
					rytmp = rys(l, k)
					rztmp = rzs(l, k)

					boolCentral = isCentral(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
					if(boolCentral) then

						m = m + 1
						jmolecule = (l - 1) / 29
						SnnCentral(m, 1, i) = k
						SnnCentral(m, 2, i) = jmolecule
						SnnCentral(m, 3, i) = l
						SnnCentral(m, 4, i) = l - 29 * jmolecule
						ictSnnCentral(i) = ictSnnCentral(i) + 1

					endif

					bool26 = isnearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
					if(bool26.and.(.not.boolCentral)) then

						n = n + 1
						jmolecule = (l - 1) / 29
						Snn26(n, 1, i) = k
						Snn26(n, 2, i) = jmolecule
						Snn26(n, 3, i) = l
						Snn26(n, 4, i) = l - 29 * jmolecule
						ictSnn26(i) = ictSnn26(i) + 1

					endif

				enddo

			enddo
			
			! north boxes (constant +y)
			ibstart = ib1000 + 100 * j + 20 - 2
			ibend = ibstart + 4

			do k = ibstart, ibend
			!if(iaq.eq.13.and.i.eq.1) write(*,*) ib1000, k

				ictk = 29 * icts(k)

				do l = 1, ictk

					rxtmp = rxs(l, k)
					rytmp = rys(l, k)
					rztmp = rzs(l, k)

					boolCentral = isCentral(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
					if(boolCentral) then

						m = m + 1
						jmolecule = (l - 1) / 29
						SnnCentral(m, 1, i) = k
						SnnCentral(m, 2, i) = jmolecule
						SnnCentral(m, 3, i) = l
						SnnCentral(m, 4, i) = l - 29 * jmolecule
						ictSnnCentral(i) = ictSnnCentral(i) + 1

					endif

					bool26 = isnearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
					if(bool26.and.(.not.boolCentral)) then

						n = n + 1
						jmolecule = (l - 1) / 29
						Snn26(n, 1, i) = k
						Snn26(n, 2, i) = jmolecule
						Snn26(n, 3, i) = l
						Snn26(n, 4, i) = l - 29 * jmolecule
						ictSnn26(i) = ictSnn26(i) + 1

					endif

				enddo

			enddo
			
			! west boxes (constant -x)
			ibstart = ib1000 + 100 * j - 10 - 2
			ibend = ibstart + 20

			do k = ibstart, ibend, 10
			!if(iaq.eq.13.and.i.eq.1) write(*,*) ib1000, k

				ictk = 29 * icts(k)

				do l = 1, ictk

					rxtmp = rxs(l, k)
					rytmp = rys(l, k)
					rztmp = rzs(l, k)
					
					boolCentral = isCentral(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
					if(boolCentral) then

						m = m + 1
						jmolecule = (l - 1) / 29
						SnnCentral(m, 1, i) = k
						SnnCentral(m, 2, i) = jmolecule
						SnnCentral(m, 3, i) = l
						SnnCentral(m, 4, i) = l - 29 * jmolecule
						ictSnnCentral(i) = ictSnnCentral(i) + 1

					endif

					bool26 = isnearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
					if(bool26.and.(.not.boolCentral)) then

						n = n + 1
						jmolecule = (l - 1) / 29
						Snn26(n, 1, i) = k
						Snn26(n, 2, i) = jmolecule
						Snn26(n, 3, i) = l
						Snn26(n, 4, i) = l - 29 * jmolecule
						ictSnn26(i) = ictSnn26(i) + 1

					endif
					
				enddo

			enddo

			! east boxes (constant +x)
			ibstart = ib1000 + 100 * j - 10 + 2
			ibend = ibstart + 20

			do k = ibstart, ibend, 10
			!if(iaq.eq.13.and.i.eq.1) write(*,*) ib1000, k

				ictk = 29 * icts(k)

				do l = 1, ictk

					rxtmp = rxs(l, k)
					rytmp = rys(l, k)
					rztmp = rzs(l, k)

					boolCentral = isCentral(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
					if(boolCentral) then

						m = m + 1
						jmolecule = (l - 1) / 29
						SnnCentral(m, 1, i) = k
						SnnCentral(m, 2, i) = jmolecule
						SnnCentral(m, 3, i) = l
						SnnCentral(m, 4, i) = l - 29 * jmolecule
						ictSnnCentral(i) = ictSnnCentral(i) + 1

					endif

					bool26 = isnearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)
					if(bool26.and.(.not.boolCentral)) then

						n = n + 1
						jmolecule = (l - 1) / 29
						Snn26(n, 1, i) = k
						Snn26(n, 2, i) = jmolecule
						Snn26(n, 3, i) = l
						Snn26(n, 4, i) = l - 29 * jmolecule
						ictSnn26(i) = ictSnn26(i) + 1

					endif

				enddo

			enddo

		enddo

		if(ictSnnCentral(i).gt.ictmaxSnn) write(*,*) "Bound exceeded in buildSnn."
		if(ictSnn26(i).gt.ictmaxSnn) write(*,*) "Bound exceeded in buildSnn."
		!if(iaq.eq.13) write(*,*) "COUNT",i,ictSnnCentral(i),ictSnn26(i)

	enddo

	return
	end subroutine
