	subroutine east_west_Onm4(comm3Dq, iaq, ibq, icq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"

	integer*4 comm3Dq, iaq, ibq, icq, j, k, l
	integer*4 ib512, ib64to512
	integer*4 itag, ierr, istatus(mpi_status_size)
	complex*16 csbuf(kewOnm4), crbuf(kewOnm4)

	itag = 0

	! send EAST receive from WEST
	l = 0
	do k = 3, 63, 4

		ib512 = ib64to512(k)

		do j = 1, lp

			l = l + 1
			csbuf(l) = cO4(j, ib512)

		enddo

	enddo
	
	do k = 4, 64, 4

		ib512 = ib64to512(k)

		do j = 1, lp

			l = l + 1
			csbuf(l) = cO4(j, ib512)

		enddo

	enddo

	call mpi_sendrecv(csbuf, kewOnm4, mpi_complex16, ibq, itag, &
	crbuf, kewOnm4, mpi_complex16, icq, itag, comm3Dq, istatus, ierr)

	if(icq.ne.mpi_proc_null) then

		l = 0
		do k = 3, 63, 4

			ib512 = ib64to512(k) - 4

			do j = 1, lp

				l = l + 1
				cO4(j, ib512) = crbuf(l)

			enddo

		enddo
		
		do k = 4, 64, 4

			ib512 = ib64to512(k) - 4

			do j = 1, lp

				l = l + 1
				cO4(j, ib512) = crbuf(l)

			enddo

		enddo

		if(l.ne.kewOnm4) write(*,*) "ewOnm4"

	endif
	
	! send WEST receive from EAST
	l = 0
	do k = 2, 62, 4

		ib512 = ib64to512(k)

		do j = 1, lp

			l = l + 1
			csbuf(l) = cO4(j, ib512)

		enddo

	enddo
	
	do k = 1, 61, 4

		ib512 = ib64to512(k)

		do j = 1, lp

			l = l + 1
			csbuf(l) = cO4(j, ib512)

		enddo

	enddo

	call mpi_sendrecv(csbuf, kewOnm4, mpi_complex16, icq, itag, &
	crbuf, kewOnm4, mpi_complex16, ibq, itag, comm3Dq, istatus, ierr)
	
	if(ibq.ne.mpi_proc_null) then

		l = 0
		do k = 2, 62, 4

			ib512 = ib64to512(k) + 4

			do j = 1, lp

				l = l + 1
				cO4(j, ib512) = crbuf(l)

			enddo

		enddo
		
		do k = 1, 61, 4

			ib512 = ib64to512(k) + 4

			do j = 1, lp

				l = l + 1
				cO4(j, ib512) = crbuf(l)

			enddo

		enddo

		if(l.ne.kewOnm4) write(*,*) "ewOnm4"

	endif

	return
	end subroutine
