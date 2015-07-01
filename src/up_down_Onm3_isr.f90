	subroutine up_down_Onm3_isr(comm3Dqq, iaqq, ibqq, icqq, idqq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"

	integer*4 comm3Dqq, iaqq, ibqq, icqq, idqq(3)
	integer*4 j, l, ib216
	integer*4 itag, ierr, ireq(4), istatus(mpi_status_size, 4)
	complex*16 csbuf(kudOnm3), crbuf(kudOnm3), crbuf2(kudOnm3)

	itag = 0

	! send UP send DOWN
	l = 0
	do ib216 = 73, 144

		do j = 1, lp

			l = l + 1
			csbuf(l) = cO3(j, ib216)

		enddo

	enddo
	
	call mpi_isend(csbuf, kudOnm3, mpi_complex16, ibqq, itag, comm3Dqq, ireq(1), ierr)
	call mpi_isend(csbuf, kudOnm3, mpi_complex16, icqq, itag, comm3Dqq, ireq(2), ierr)
	call mpi_irecv(crbuf, kudOnm3, mpi_complex16, icqq, itag, comm3Dqq, ireq(3), ierr)
	call mpi_irecv(crbuf2, kudOnm3, mpi_complex16, ibqq, itag, comm3Dqq, ireq(4), ierr)
	
	call computeOnm2(iaqq)
	call ewnsud_Onm2(comm3Dqq)
	call computeLjk2(iaqq, idqq)
	call shiftLjk2(iaqq)
	
	call mpi_waitall(4, ireq, istatus, ierr)
	
	! receive from DOWN
	if(icqq.ne.mpi_proc_null) then

		l = 0
		do ib216 = 1, 72

			do j = 1, lp

				l = l + 1
				cO3(j, ib216) = crbuf(l)

			enddo

		enddo

	endif

	! receive from UP
	if(ibqq.ne.mpi_proc_null) then

		l = 0
		do ib216 = 145, 216

			do j = 1, lp

				l = l + 1
				cO3(j, ib216) = crbuf2(l)

			enddo

		enddo

	endif

	return
	end subroutine
