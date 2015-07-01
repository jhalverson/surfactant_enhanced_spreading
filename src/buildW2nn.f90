	! could just loop over H atoms, no need to store type since it is H
	! check 5^3 - 3^3 = 98 boxes

	subroutine buildW2nn(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 i, j, k, l, m, ib512, ib64to512
	integer*4 ibstart, ibend, ictk, jmolecule

	real*8 rxtmp, rytmp, rztmp
	real*8 rctr4x, rctr4y, rctr4z

	logical bool26, isNearest

	do i = 1, 64

		ictW2nn(i) = 0
		m = 0

		ib512 = ib64to512(i)
		rctr4x = rcen4x(i)
		rctr4y = rcen4y(i)
		rctr4z = rcen4z(i)

		! lower 5 by 5
		do j = 1, 5

			ibstart = ib512 - 128 - 18 + (j - 1) * 8
			ibend = ibstart + 4

			do k = ibstart, ibend

				ictk = 3 * ictw2(k)

				do l = 1, ictk

					rxtmp = rxw2(l, k)
					rytmp = ryw2(l, k)
					rztmp = rzw2(l, k)

					bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)

					if(bool26) then

						m = m + 1
						jmolecule = (l - 1) / 3
						W2nn(m, 1, i) = k
						W2nn(m, 3, i) = l
						W2nn(m, 4, i) = l - 3 * jmolecule
						ictW2nn(i) = ictW2nn(i) + 1

					endif

				enddo

			enddo

		enddo
		
		! upper 5 by 5
		do j = 1, 5

			ibstart = ib512 + 128 - 18 + (j - 1) * 8
			ibend = ibstart + 4

			do k = ibstart, ibend

				ictk = 3 * ictw2(k)

				do l = 1, ictk

					rxtmp = rxw2(l, k)
					rytmp = ryw2(l, k)
					rztmp = rzw2(l, k)

					bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)

					if(bool26) then

						m = m + 1
						jmolecule = (l - 1) / 3
						W2nn(m, 1, i) = k
						W2nn(m, 3, i) = l
						W2nn(m, 4, i) = l - 3 * jmolecule
						ictW2nn(i) = ictW2nn(i) + 1

					endif

				enddo

			enddo

		enddo
		
		! square rings
		do j = -1, 1

			! south boxes (constant -y)
			ibstart = ib512 + 64 * j - 16 - 2
			ibend = ibstart + 4

			do k = ibstart, ibend

				ictk = 3 * ictw2(k)

				do l = 1, ictk

					rxtmp = rxw2(l, k)
					rytmp = ryw2(l, k)
					rztmp = rzw2(l, k)

					bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)

					if(bool26) then

						m = m + 1
						jmolecule = (l - 1) / 3
						W2nn(m, 1, i) = k
						W2nn(m, 3, i) = l
						W2nn(m, 4, i) = l - 3 * jmolecule
						ictW2nn(i) = ictW2nn(i) + 1

					endif

				enddo

			enddo

			! north boxes (constant +y)
			ibstart = ib512 + 64 * j + 16 - 2
			ibend = ibstart + 4

			do k = ibstart, ibend

				ictk = 3 * ictw2(k)

				do l = 1, ictk

					rxtmp = rxw2(l, k)
					rytmp = ryw2(l, k)
					rztmp = rzw2(l, k)

					bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)

					if(bool26) then

						m = m + 1
						jmolecule = (l - 1) / 3
						W2nn(m, 1, i) = k
						W2nn(m, 3, i) = l
						W2nn(m, 4, i) = l - 3 * jmolecule
						ictW2nn(i) = ictW2nn(i) + 1

					endif

				enddo

			enddo

			! west boxes (constant -x)
			ibstart = ib512 + 64 * j - 8 - 2
			ibend = ibstart + 16

			do k = ibstart, ibend, 8

				ictk = 3 * ictw2(k)

				do l = 1, ictk

					rxtmp = rxw2(l, k)
					rytmp = ryw2(l, k)
					rztmp = rzw2(l, k)

					bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)

					if(bool26) then

						m = m + 1
						jmolecule = (l - 1) / 3
						W2nn(m, 1, i) = k
						W2nn(m, 3, i) = l
						W2nn(m, 4, i) = l - 3 * jmolecule
						ictW2nn(i) = ictW2nn(i) + 1

					endif

				enddo

			enddo

			! east boxes (constant +x)
			ibstart = ib512 + 64 * j - 8 + 2
			ibend = ibstart + 16

			do k = ibstart, ibend, 8

				ictk = 3 * ictw2(k)

				do l = 1, ictk

					rxtmp = rxw2(l, k)
					rytmp = ryw2(l, k)
					rztmp = rzw2(l, k)

					bool26 = isNearest(rxtmp, rytmp, rztmp, rctr4x, rctr4y, rctr4z)

					if(bool26) then

						m = m + 1
						jmolecule = (l - 1) / 3
						W2nn(m, 1, i) = k
						W2nn(m, 3, i) = l
						W2nn(m, 4, i) = l - 3 * jmolecule
						ictW2nn(i) = ictW2nn(i) + 1

					endif

				enddo

			enddo

		enddo

		if(ictW2nn(i).gt.ictmaxW2nn) write(*,*) "buildW2nn ictmaxW2nn"

	enddo

	return
	end subroutine
