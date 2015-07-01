	subroutine north_south_Onm3(comm3Dq, iaq, ibq, icq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"

	integer*4 comm3Dq, iaq, ibq, icq
	integer*4 j, k, l, ib216
	integer*4 itag, ierr, istatus(mpi_status_size)
	complex*16 csbuf(knsOnm3), crbuf(knsOnm3)

	itag = 0

	l = 0
	do ib216 = 85, 96

		do j = 1, lp

			l = l + 1
			csbuf(l) = cO3(j, ib216)

		enddo

	enddo

	do ib216 = 121, 132

		do j = 1, lp

			l = l + 1
			csbuf(l) = cO3(j, ib216)

		enddo

	enddo

	! send NORTH receive from SOUTH
	call mpi_sendrecv(csbuf, knsOnm3, mpi_complex16, ibq, itag, &
	crbuf, knsOnm3, mpi_complex16, icq, itag, comm3Dq, istatus, ierr)

	if(icq.ne.mpi_proc_null) then

		l = 0
		do ib216 = 73, 84

			do j = 1, lp

				l = l + 1
				cO3(j, ib216) = crbuf(l)

			enddo

		enddo

		do ib216 = 109, 120

			do j = 1, lp

				l = l + 1
				cO3(j, ib216) = crbuf(l)

			enddo

		enddo

	endif

	! send SOUTH receive from NORTH
	call mpi_sendrecv(csbuf, knsOnm3, mpi_complex16, icq, itag, &
	crbuf, knsOnm3, mpi_complex16, ibq, itag, comm3Dq, istatus, ierr)

	if(ibq.ne.mpi_proc_null) then

		l = 0
		do ib216 = 97, 108

			do j = 1, lp

				l = l + 1
				cO3(j, ib216) = crbuf(l)

			enddo

		enddo

		do ib216 = 133, 144

			do j = 1, lp

				l = l + 1
				cO3(j, ib216) = crbuf(l)

			enddo

		enddo

	endif

	return
	end subroutine
