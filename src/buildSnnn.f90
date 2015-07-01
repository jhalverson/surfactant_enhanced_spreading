	! check 7^3 - 5^3 = 218 boxes
	subroutine buildSnnn(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 i, j, k, l, m, ib1000, ib64to1000
	integer*4 ibstart, ibend, ictk, jmolecule

	real*8 rxtmp, rytmp, rztmp
	real*8 rctr4x, rctr4y, rctr4z

	logical bool26, isNearest

	do i = 1, 64

		ictSnnn(i) = 0
		m = 0

		ib1000 = ib64to1000(i)
		rctr4x = rcen4x(i)
		rctr4y = rcen4y(i)
		rctr4z = rcen4z(i)

		! lower 7 by 7
		do j = 1, 7

			ibstart = ib1000 - 300 - 33 + (j - 1) * 10
			ibend = ibstart + 6

			do k = ibstart, ibend

				ictk = 29 * icts(k)

				do l = 1, ictk

					rxtmp = rxs(l, k)
					rytmp = rys(l, k)
					rztmp = rzs(l, k)

					bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)

					if(bool26) then

						m = m + 1
						jmolecule = (l - 1) / 29
						Snnn(m, 1, i) = k
						Snnn(m, 3, i) = l
						Snnn(m, 4, i) = l - 29 * jmolecule
						ictSnnn(i) = ictSnnn(i) + 1

					endif

				enddo

			enddo

		enddo

		! upper 7 by 7
		do j = 1, 7

			ibstart = ib1000 + 300 - 33 + (j - 1) * 10
			ibend = ibstart + 6

			do k = ibstart, ibend

				ictk = 29 * icts(k)

				do l = 1, ictk

					rxtmp = rxs(l, k)
					rytmp = rys(l, k)
					rztmp = rzs(l, k)

					bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)

					if(bool26) then

						m = m + 1
						jmolecule = (l - 1) / 29
						Snnn(m, 1, i) = k
						Snnn(m, 3, i) = l
						Snnn(m, 4, i) = l - 29 * jmolecule
						ictSnnn(i) = ictSnnn(i) + 1

					endif

				enddo

			enddo

		enddo

		! 5 square rings
		do j = -2, 2

			! south 7 boxes (constant -y)
			ibstart = ib1000 + 100 * j - 30 - 3
			ibend = ibstart + 6

			do k = ibstart, ibend

				ictk = 29 * icts(k)

				do l = 1, ictk

					rxtmp = rxs(l, k)
					rytmp = rys(l, k)
					rztmp = rzs(l, k)

					bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)

					if(bool26) then

						m = m + 1
						jmolecule = (l - 1) / 29
						Snnn(m, 1, i) = k
						Snnn(m, 3, i) = l
						Snnn(m, 4, i) = l - 29 * jmolecule
						ictSnnn(i) = ictSnnn(i) + 1

					endif

				enddo

			enddo

			! north 7 boxes (constant +y)
			ibstart = ib1000 + 100 * j + 30 - 3
			ibend = ibstart + 6

			do k = ibstart, ibend

				ictk = 29 * icts(k)

				do l = 1, ictk

					rxtmp = rxs(l, k)
					rytmp = rys(l, k)
					rztmp = rzs(l, k)

					bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)

					if(bool26) then

						m = m + 1
						jmolecule = (l - 1) / 29
						Snnn(m, 1, i) = k
						Snnn(m, 3, i) = l
						Snnn(m, 4, i) = l - 29 * jmolecule
						ictSnnn(i) = ictSnnn(i) + 1

					endif

				enddo

			enddo

			! west 7 boxes (constant -x)
			ibstart = ib1000 + 100 * j - 20 - 3
			ibend = ibstart + 40

			do k = ibstart, ibend, 10

				ictk = 29 * icts(k)

				do l = 1, ictk

					rxtmp = rxs(l, k)
					rytmp = rys(l, k)
					rztmp = rzs(l, k)

					bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)

					if(bool26) then

						m = m + 1
						jmolecule = (l - 1) / 29
						Snnn(m, 1, i) = k
						Snnn(m, 3, i) = l
						Snnn(m, 4, i) = l - 29 * jmolecule
						ictSnnn(i) = ictSnnn(i) + 1

					endif

				enddo

			enddo

			! east 7 boxes (constant +x)
			ibstart = ib1000 + 100 * j - 20 + 3
			ibend = ibstart + 40

			do k = ibstart, ibend, 10

				ictk = 29 * icts(k)

				do l = 1, ictk

					rxtmp = rxs(l, k)
					rytmp = rys(l, k)
					rztmp = rzs(l, k)

					bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)

					if(bool26) then

						m = m + 1
						jmolecule = (l - 1) / 29
						Snnn(m, 1, i) = k
						Snnn(m, 3, i) = l
						Snnn(m, 4, i) = l - 29 * jmolecule
						ictSnnn(i) = ictSnnn(i) + 1

					endif

				enddo

			enddo

		enddo

		if(ictSnnn(i).gt.ictmaxSnnn) write(*,*) "buildSnnn bound"

	enddo

	return
	end subroutine
