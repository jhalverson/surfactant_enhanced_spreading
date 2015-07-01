	subroutine northSR(comm3Dq, iaq, ibq, icq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"

	integer*4 comm3Dq, iaq, ibq, icq
	integer*4 j, k, l, m, n
	integer*4 ib1000, ib512, ib216, ictmolecule
	integer*4 i3l, ictcrd
	integer*4 iloop, ibstart, ibend
	integer*4 itag, ierr, istatus(mpi_status_size)

	integer*4 ictbuf(208), ictrbuf(208)
	real*8 rcrdbuf(knsRws), rcrdrbuf(knsRws)

	! load arrays for sending

	k = 0
	l = 0
	do iloop = 1, 4

		ibstart = 341 + (iloop - 1) * 100
		ibend = 350 + (iloop - 1) * 100

		do ib1000 = ibstart, ibend

			k = k + 1
			ictbuf(k) = icts(ib1000)
			ictmolecule = icts(ib1000)

			do j = 1, ictmolecule

				do m = 1, 29

					n = (j - 1) * 29 + m
					rcrdbuf(l + 1) = rxs(n, ib1000)
					rcrdbuf(l + 2) = rys(n, ib1000)
					rcrdbuf(l + 3) = rzs(n, ib1000)
					l = l + 3

				enddo

			enddo

		enddo

	enddo

	do iloop = 1, 4

		ibstart = 351 + (iloop - 1) * 100
		ibend = 360 + (iloop - 1) * 100

		do ib1000 = ibstart, ibend

			k = k + 1
			ictbuf(k) = icts(ib1000)
			ictmolecule = icts(ib1000)

			do j = 1, ictmolecule

				do m = 1, 29

					n = (j - 1) * 29 + m
					rcrdbuf(l + 1) = rxs(n, ib1000)
					rcrdbuf(l + 2) = rys(n, ib1000)
					rcrdbuf(l + 3) = rzs(n, ib1000)
					l = l + 3

				enddo

			enddo

		enddo

	enddo

	do iloop = 1, 4

		ibstart = 161 + (iloop - 1) * 64
		ibend = 168 + (iloop - 1) * 64

		do ib512 = ibstart, ibend

			k = k + 1
			ictbuf(k) = ictw2(ib512)
			ictmolecule = ictw2(ib512)

			do j = 1, ictmolecule

				do m = 1, 3

					n = (j - 1) * 3 + m
					rcrdbuf(l + 1) = rxw2(n, ib512)
					rcrdbuf(l + 2) = ryw2(n, ib512)
					rcrdbuf(l + 3) = rzw2(n, ib512)
					l = l + 3

				enddo

			enddo

		enddo

	enddo

	do iloop = 1, 4

		ibstart = 361 + (iloop - 1) * 100
		ibend = 370 + (iloop - 1) * 100

		do ib1000 = ibstart, ibend

			k = k + 1
			ictbuf(k) = icts(ib1000)
			ictmolecule = icts(ib1000)

			do j = 1, ictmolecule

				do m = 1, 29

					n = (j - 1) * 29 + m
					rcrdbuf(l + 1) = rxs(n, ib1000)
					rcrdbuf(l + 2) = rys(n, ib1000)
					rcrdbuf(l + 3) = rzs(n, ib1000)
					l = l + 3

				enddo

			enddo

		enddo

	enddo

	do iloop = 1, 4

		ibstart = 169 + (iloop - 1) * 64
		ibend = 176 + (iloop - 1) * 64

		do ib512 = ibstart, ibend

			k = k + 1
			ictbuf(k) = ictw2(ib512)
			ictmolecule = ictw2(ib512)

			do j = 1, ictmolecule

				do m = 1, 3

					n = (j - 1) * 3 + m
					rcrdbuf(l + 1) = rxw2(n, ib512)
					rcrdbuf(l + 2) = ryw2(n, ib512)
					rcrdbuf(l + 3) = rzw2(n, ib512)
					l = l + 3

				enddo

			enddo

		enddo

	enddo

	do iloop = 1, 4

		ibstart = 61 + (iloop - 1) * 36
		ibend = 66 + (iloop - 1) * 36

		do ib216 = ibstart, ibend

			k = k + 1
			ictbuf(k) = ictw1(ib216)
			ictmolecule = ictw1(ib216)

			do j = 1, ictmolecule

				do m = 1, 3

					n = (j - 1) * 3 + m
					rcrdbuf(l + 1) = rxw1(n, ib216)
					rcrdbuf(l + 2) = ryw1(n, ib216)
					rcrdbuf(l + 3) = rzw1(n, ib216)
					l = l + 3

				enddo

			enddo

		enddo

	enddo

	! check
	if(l.gt.knsRws) write(*,*) "northSR knsRws"
	if(k.ne.208) write(*,*) "northSR k load"

	itag = 0
	i3l = l

	! send north receive south
	
	call mpi_sendrecv(ictbuf, 208, mpi_integer4, ibq, i3l, ictrbuf, 208, mpi_integer4, &
	icq, mpi_any_tag, comm3Dq, istatus, ierr)
	
	ictcrd = istatus(mpi_tag)
	if(ictcrd.lt.0) ictcrd = 0
	
	call mpi_sendrecv(rcrdbuf, i3l, mpi_real8, ibq, itag, rcrdrbuf, ictcrd, mpi_real8, &
	icq, itag, comm3Dq, istatus, ierr)
	
	! unload arrays if received
	
	if(icq.ne.mpi_proc_null) then

		k = 0
		l = 0
		do iloop = 1, 4

			ibstart = 301 + (iloop - 1) * 100
			ibend = 310 + (iloop - 1) * 100

			do ib1000 = ibstart, ibend

				k = k + 1
				icts(ib1000) = ictrbuf(k)
				ictmolecule = icts(ib1000)

				do j = 1, ictmolecule

					do m = 1, 29

						n = (j - 1) * 29 + m
						rxs(n, ib1000) = rcrdrbuf(l + 1)
						rys(n, ib1000) = rcrdrbuf(l + 2)
						rzs(n, ib1000) = rcrdrbuf(l + 3)
						l = l + 3

					enddo

				enddo

			enddo

		enddo

		do iloop = 1, 4

			ibstart = 311 + (iloop - 1) * 100
			ibend = 320 + (iloop - 1) * 100

			do ib1000 = ibstart, ibend

				k = k + 1
				icts(ib1000) = ictrbuf(k)
				ictmolecule = icts(ib1000)

				do j = 1, ictmolecule

					do m = 1, 29

						n = (j - 1) * 29 + m
						rxs(n, ib1000) = rcrdrbuf(l + 1)
						rys(n, ib1000) = rcrdrbuf(l + 2)
						rzs(n, ib1000) = rcrdrbuf(l + 3)
						l = l + 3

					enddo

				enddo

			enddo

		enddo

		do iloop = 1, 4

			ibstart = 129 + (iloop - 1) * 64
			ibend = 136 + (iloop - 1) * 64

			do ib512 = ibstart, ibend

				k = k + 1
				ictw2(ib512) = ictrbuf(k)
				ictmolecule = ictw2(ib512)

				do j = 1, ictmolecule

					do m = 1, 3

						n = (j - 1) * 3 + m
						rxw2(n, ib512) = rcrdrbuf(l + 1)
						ryw2(n, ib512) = rcrdrbuf(l + 2)
						rzw2(n, ib512) = rcrdrbuf(l + 3)
						l = l + 3

					enddo

				enddo

			enddo

		enddo

		do iloop = 1, 4

			ibstart = 321 + (iloop - 1) * 100
			ibend = 330 + (iloop - 1) * 100

			do ib1000 = ibstart, ibend

				k = k + 1
				icts(ib1000) = ictrbuf(k)
				ictmolecule = icts(ib1000)

				do j = 1, ictmolecule

					do m = 1, 29

						n = (j - 1) * 29 + m
						rxs(n, ib1000) = rcrdrbuf(l + 1)
						rys(n, ib1000) = rcrdrbuf(l + 2)
						rzs(n, ib1000) = rcrdrbuf(l + 3)
						l = l + 3

					enddo

				enddo

			enddo

		enddo

		do iloop = 1, 4

			ibstart = 137 + (iloop - 1) * 64
			ibend = 144 + (iloop - 1) * 64

			do ib512 = ibstart, ibend

				k = k + 1
				ictw2(ib512) = ictrbuf(k)
				ictmolecule = ictw2(ib512)

				do j = 1, ictmolecule

					do m = 1, 3

						n = (j - 1) * 3 + m
						rxw2(n, ib512) = rcrdrbuf(l + 1)
						ryw2(n, ib512) = rcrdrbuf(l + 2)
						rzw2(n, ib512) = rcrdrbuf(l + 3)
						l = l + 3

					enddo

				enddo

			enddo

		enddo

		do iloop = 1, 4

			ibstart = 37 + (iloop - 1) * 36
			ibend = 42 + (iloop - 1) * 36

			do ib216 = ibstart, ibend

				k = k + 1
				ictw1(ib216) = ictrbuf(k)
				ictmolecule = ictw1(ib216)

				do j = 1, ictmolecule

					do m = 1, 3

						n = (j - 1) * 3 + m
						rxw1(n, ib216) = rcrdrbuf(l + 1)
						ryw1(n, ib216) = rcrdrbuf(l + 2)
						rzw1(n, ib216) = rcrdrbuf(l + 3)
						l = l + 3

					enddo

				enddo

			enddo

		enddo

	if(l.ne.ictcrd) write(*,*) "northSR l"
	if(k.ne.208) write(*,*) "northSR k unload"

	endif

	return
	end subroutine
