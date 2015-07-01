	subroutine north_south_Onm4(comm3Dq, iaq, ibq, icq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"

	integer*4 comm3Dq, iaq, ibq, icq
	integer*4 ib512, ibstart, ibend, iloop, j, k, l
	integer*4 itag, ierr, istatus(mpi_status_size)
	complex*16 csbuf(knsOnm4), crbuf(knsOnm4)

	itag = 0

	! send NORTH receive from SOUTH
	l = 0
	do iloop = 1, 4

		ibstart = 161 + (iloop - 1) * 64
		ibend = 168 + (iloop - 1) * 64

		do ib512 = ibstart, ibend

			do j = 1, lp

				l = l + 1
				csbuf(l) = cO4(j, ib512)

			enddo

		enddo

	enddo

	do iloop = 1, 4

		ibstart = 169 + (iloop - 1) * 64
		ibend = 176 + (iloop - 1) * 64

		do ib512 = ibstart, ibend

			do j = 1, lp

				l = l + 1
				csbuf(l) = cO4(j, ib512)

			enddo

		enddo

	enddo

	call mpi_sendrecv(csbuf, knsOnm4, mpi_complex16, ibq, itag, &
	crbuf, knsOnm4, mpi_complex16, icq, itag, comm3Dq, istatus, ierr)

	if(icq.ne.mpi_proc_null) then

		l = 0
		do iloop = 1, 4

			ibstart = 129 + (iloop - 1) * 64
			ibend = 136 + (iloop - 1) * 64
			
			do ib512 = ibstart, ibend
			
				do j = 1, lp
				
					l = l + 1
					cO4(j, ib512) = crbuf(l)
					
				enddo
				
			enddo
			
		enddo
		
		do iloop = 1, 4

			ibstart = 137 + (iloop - 1) * 64
			ibend = 144 + (iloop - 1) * 64
			
			do ib512 = ibstart, ibend
			
				do j = 1, lp
				
					l = l + 1
					cO4(j, ib512) = crbuf(l)
					
				enddo
				
			enddo
			
		enddo

		if(l.ne.knsOnm4) write(*,*) "nsOnm4"
		
	endif
	
	! send SOUTH receive from NORTH
	l = 0
	do iloop = 1, 4

		ibstart = 153 + (iloop - 1) * 64
		ibend = 160 + (iloop - 1) * 64

		do ib512 = ibstart, ibend

			do j = 1, lp

				l = l + 1
				csbuf(l) = cO4(j, ib512)

			enddo

		enddo

	enddo

	do iloop = 1, 4

		ibstart = 145 + (iloop - 1) * 64
		ibend = 152 + (iloop - 1) * 64

		do ib512 = ibstart, ibend

			do j = 1, lp

				l = l + 1
				csbuf(l) = cO4(j, ib512)

			enddo

		enddo

	enddo

	call mpi_sendrecv(csbuf, knsOnm4, mpi_complex16, icq, itag, &
	crbuf, knsOnm4, mpi_complex16, ibq, itag, comm3Dq, istatus, ierr)

	if(ibq.ne.mpi_proc_null) then

		l = 0
		do iloop = 1, 4

			ibstart = 185 + (iloop - 1) * 64
			ibend = 192 + (iloop - 1) * 64
			
			do ib512 = ibstart, ibend
			
				do j = 1, lp
				
					l = l + 1
					cO4(j, ib512) = crbuf(l)
					
				enddo
				
			enddo
			
		enddo
		
		do iloop = 1, 4

			ibstart = 177 + (iloop - 1) * 64
			ibend = 184 + (iloop - 1) * 64
			
			do ib512 = ibstart, ibend
			
				do j = 1, lp
				
					l = l + 1
					cO4(j, ib512) = crbuf(l)
					
				enddo
				
			enddo
			
		enddo

		if(l.ne.knsOnm4) write(*,*) "nsOnm4"
		
	endif

	return
	end subroutine
