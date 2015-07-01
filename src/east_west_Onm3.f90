	subroutine east_west_Onm3(comm3Dq, iaq, ibq, icq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"

	integer*4 comm3Dq, iaq, ibq, icq
	integer*4 j, k, l, ib216, ib8to216
	integer*4 itag, ierr, istatus(mpi_status_size)
	complex*16 csbuf(kewOnm3), crbuf(kewOnm3)

	itag = 0

	l = 0
	do k = 1, 8

		ib216 = ib8to216(k)

		do j = 1, lp

			l = l + 1
			csbuf(l) = cO3(j, ib216)

		enddo

	enddo

	! send EAST receive from WEST
	call mpi_sendrecv(csbuf, kewOnm3, mpi_complex16, ibq, itag, &
	crbuf, kewOnm3, mpi_complex16, icq, itag, comm3Dq, istatus, ierr)

	if(icq.ne.mpi_proc_null) then

		l = 0
		do k = 1, 8

			ib216 = ib8to216(k) - 2

			do j = 1, lp

				l = l + 1
				cO3(j, ib216) = crbuf(l)

			enddo

		enddo

	endif
	
	! send WEST receive from EAST
	call mpi_sendrecv(csbuf, kewOnm3, mpi_complex16, icq, itag, &
	crbuf, kewOnm3, mpi_complex16, ibq, itag, comm3Dq, istatus, ierr)

	if(ibq.ne.mpi_proc_null) then

		l = 0
		do k = 1, 8

			ib216 = ib8to216(k) + 2

			do j = 1, lp

				l = l + 1
				cO3(j, ib216) = crbuf(l)

			enddo

		enddo

	endif

	return
	end subroutine
