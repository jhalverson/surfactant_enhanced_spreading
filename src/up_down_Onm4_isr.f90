	subroutine up_down_Onm4_isr(comm3Dq, iaq, ibq, icq, idq, ieq, ifq, igq, ihq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"
	
	integer*4 comm3Dq, iaq, ibq, icq, idq, ieq, ifq, igq, ihq(3)
	integer*4 ib512, j, l
	integer*4 itag, ierr, ireq(4), istatus(mpi_status_size, 4)
	complex*16 csbuf(kudOnm4), csbuf2(kudOnm4)
	complex*16 crbuf(kudOnm4), crbuf2(kudOnm4)

	itag = 0

	! send UP
	l = 0
	do ib512 = 257, 384

		do j = 1, lp

			l = l + 1
			csbuf(l) = cO4(j, ib512)

		enddo

	enddo
	
	! send DOWN
	l = 0
	do ib512 = 129, 256

		do j = 1, lp

			l = l + 1
			csbuf2(l) = cO4(j, ib512)

		enddo

	enddo

	call mpi_isend(csbuf, kudOnm4, mpi_complex16, ifq, itag, comm3Dq, ireq(1), ierr)
	call mpi_isend(csbuf2, kudOnm4, mpi_complex16, igq, itag, comm3Dq, ireq(2), ierr)
	call mpi_irecv(crbuf, kudOnm4, mpi_complex16, igq, itag, comm3Dq, ireq(3), ierr)
	call mpi_irecv(crbuf2, kudOnm4, mpi_complex16, ifq, itag, comm3Dq, ireq(4), ierr)

	call computeOnm3(iaq)
	call east_west_Onm3(comm3Dq, iaq, ibq, icq)
	call north_south_Onm3(comm3Dq, iaq, idq, ieq)
	call up_down_Onm3_isr(comm3Dq, iaq, ifq, igq, ihq)
	call computeLjk3(iaq)
	call shiftLjk3(iaq)

	call mpi_waitall(4, ireq, istatus, ierr)
	
	! receive from DOWN
	if(igq.ne.mpi_proc_null) then

		l = 0
		do ib512 = 1, 128

			do j = 1, lp

				l = l + 1
				cO4(j, ib512) = crbuf(l)

			enddo

		enddo

	endif

	! receive from UP
	if(ifq.ne.mpi_proc_null) then

		l = 0
		do ib512 = 385, 512

			do j = 1, lp

				l = l + 1
				cO4(j, ib512) = crbuf2(l)

			enddo

		enddo

	endif

	return
	end subroutine
