	subroutine westSR(comm3Dq, iaq, ibq, icq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"
	
	integer*4 comm3Dq, iaq, ibq, icq
	integer*4 i, j, k, l, m, n
	integer*4 ib1000, ib64to1000, ictmolecule
	integer*4 ib512, ib64to512, ib216, ib64to216
	integer*4 i3l, ictcrd
	integer*4 itag, ierr, istatus(mpi_status_size)

	integer*4 ictbuf(96), ictrbuf(96)
	real*8 rcrdbuf(kewRws), rcrdrbuf(kewRws)

	! load arrays for sending

	k = 0
	l = 0
	do i = 3, 63, 4

		k = k + 1
		ib1000 = ib64to1000(i)
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

	do i = 2, 62, 4

		k = k + 1
		ib1000 = ib64to1000(i)
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

	do i = 2, 62, 4

		k = k + 1
		ib512 = ib64to512(i)
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

	do i = 1, 61, 4

		k = k + 1
		ib1000 = ib64to1000(i)
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

	do i = 1, 61, 4

		k = k + 1
		ib512 = ib64to512(i)
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
	
	do i = 1, 61, 4

		k = k + 1
		ib216 = ib64to216(i)
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

	! check
	if(l.gt.kewRws) write(*,*) "eastSR kewRws"
	if(k.ne.96) write(*,*) "eastSR k load"
	
	itag = 0
	i3l = l

	! send west receive from east
	
	call mpi_sendrecv(ictbuf, 96, mpi_integer4, ibq, i3l, ictrbuf, 96, mpi_integer4, &
	icq, mpi_any_tag, comm3Dq, istatus, ierr)
	
	ictcrd = istatus(mpi_tag)
	if(ictcrd.lt.0) ictcrd = 0
	
	call mpi_sendrecv(rcrdbuf, i3l, mpi_real8, ibq, itag, rcrdrbuf, ictcrd, mpi_real8, &
	icq, itag, comm3Dq, istatus, ierr)
	
	! unload arrays if received
	
	if(icq.ne.mpi_proc_null) then
		
		k = 0
		l = 0
		do i = 3, 63, 4

			k = k + 1
			ib1000 = ib64to1000(i) + 4
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

		do i = 2, 62, 4

			k = k + 1
			ib1000 = ib64to1000(i) + 4
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

		do i = 2, 62, 4

			k = k + 1
			ib512 = ib64to512(i) + 4
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

		do i = 1, 61, 4

			k = k + 1
			ib1000 = ib64to1000(i) + 4
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
		
		do i = 1, 61, 4

			k = k + 1
			ib512 = ib64to512(i) + 4
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

		do i = 1, 61, 4

			k = k + 1
			ib216 = ib64to216(i) + 4
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

		if(l.ne.ictcrd) write(*,*) "westSR l"
		if(k.ne.96) write(*,*) "westSR k unload"

	endif

	return
	end subroutine
