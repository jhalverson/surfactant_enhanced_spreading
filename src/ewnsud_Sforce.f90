	subroutine ewnsud_Sforce(comm3Dq, iaq, ibq, icq, idq, ieq, ifq, igq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"

	integer*4 comm3Dq, iaq, ibq, icq, idq, ieq, ifq, igq
	integer*4 itag, istatus(mpi_status_size), ierr
	integer*4 jstart, jend, j, iend
	integer*4 kpacket1, kpacket3, kpacket9
	integer*4 irbuf1(ksplit), irbuf3(ksplit3), irbuf9(ksplit9)
	real*8 rrbuf1(ksplit), rrbuf3(ksplit3), rrbuf9(ksplit9)

	itag = 0

	jstart = 1
	jend = jsplits
	jsplits_total = jend

	kpacket1 = jsplits

	iend = 0
	! send east receive from west
	call mpi_sendrecv(kpacket1, 1, mpi_integer4, ibq, itag, iend, 1, mpi_integer4, &
	icq, itag, comm3Dq, istatus, ierr)
	call mpi_sendrecv(isforce, kpacket1, mpi_integer4, ibq, itag, irbuf1, iend, mpi_integer4, &
	icq, itag, comm3Dq, istatus, ierr)
	call mpi_sendrecv(rsforce, kpacket1, mpi_real8, ibq, itag, rrbuf1, iend, mpi_real8, &
	icq, itag, comm3Dq, istatus, ierr)
	
	if(icq.ne.mpi_proc_null.and.iend.ne.0) then

		jstart = jend + 1
		jend = jend + iend
		jsplits_total = jend

		do j = jstart, jend

			isforce(j) = irbuf1(j - jstart + 1)
			rsforce(j) = rrbuf1(j - jstart + 1)

		enddo
	
	endif

	iend = 0
	! send west receive from east
	call mpi_sendrecv(kpacket1, 1, mpi_integer4, icq, itag, iend, 1, mpi_integer4, &
	ibq, itag, comm3Dq, istatus, ierr)
	call mpi_sendrecv(isforce, kpacket1, mpi_integer4, icq, itag, irbuf1, iend, mpi_integer4, &
	ibq, itag, comm3Dq, istatus, ierr)
	call mpi_sendrecv(rsforce, kpacket1, mpi_real8, icq, itag, rrbuf1, iend, mpi_real8, &
	ibq, itag, comm3Dq, istatus, ierr)
	
	if(ibq.ne.mpi_proc_null.and.iend.ne.0) then

		jstart = jend + 1
		jend = jend + iend
		jsplits_total = jend

		do j = jstart, jend

			isforce(j) = irbuf1(j - jstart + 1)
			rsforce(j) = rrbuf1(j - jstart + 1)

		enddo
	
	endif

	kpacket3 = jend

	iend = 0
	! send north receive from south
	call mpi_sendrecv(kpacket3, 1, mpi_integer4, idq, itag, iend, 1, mpi_integer4, &
	ieq, itag, comm3Dq, istatus, ierr)
	call mpi_sendrecv(isforce, kpacket3, mpi_integer4, idq, itag, irbuf3, iend, mpi_integer4, &
	ieq, itag, comm3Dq, istatus, ierr)
	call mpi_sendrecv(rsforce, kpacket3, mpi_real8, idq, itag, rrbuf3, iend, mpi_real8, &
	ieq, itag, comm3Dq, istatus, ierr)
	
	if(ieq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jsplits_total = jend
		
		do j = jstart, jend
		
			isforce(j) = irbuf3(j - jstart + 1)
			rsforce(j) = rrbuf3(j - jstart + 1)
		
		enddo
	
	endif
	
	iend = 0
	! send south receive from north
	call mpi_sendrecv(kpacket3, 1, mpi_integer4, ieq, itag, iend, 1, mpi_integer4, &
	idq, itag, comm3Dq, istatus, ierr)
	call mpi_sendrecv(isforce, kpacket3, mpi_integer4, ieq, itag, irbuf3, iend, mpi_integer4, &
	idq, itag, comm3Dq, istatus, ierr)
	call mpi_sendrecv(rsforce, kpacket3, mpi_real8, ieq, itag, rrbuf3, iend, mpi_real8, &
	idq, itag, comm3Dq, istatus, ierr)
	
	if(idq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jsplits_total = jend
		
		do j = jstart, jend
		
			isforce(j) = irbuf3(j - jstart + 1)
			rsforce(j) = rrbuf3(j - jstart + 1)
		
		enddo
	
	endif
	
	kpacket9 = jend
	
	iend = 0
	! send up receive from down
	call mpi_sendrecv(kpacket9, 1, mpi_integer4, ifq, itag, iend, 1, mpi_integer4, &
	igq, itag, comm3Dq, istatus, ierr)
	call mpi_sendrecv(isforce, kpacket9, mpi_integer4, ifq, itag, irbuf9, iend, mpi_integer4, &
	igq, itag, comm3Dq, istatus, ierr)
	call mpi_sendrecv(rsforce, kpacket9, mpi_real8, ifq, itag, rrbuf9, iend, mpi_real8, &
	igq, itag, comm3Dq, istatus, ierr)
	
	if(igq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jsplits_total = jend
		
		do j = jstart, jend
		
			isforce(j) = irbuf9(j - jstart + 1)
			rsforce(j) = rrbuf9(j - jstart + 1)
		
		enddo
	
	endif
	
	iend = 0
	! send down receive from up
	call mpi_sendrecv(kpacket9, 1, mpi_integer4, igq, itag, iend, 1, mpi_integer4, &
	ifq, itag, comm3Dq, istatus, ierr)
	call mpi_sendrecv(isforce, kpacket9, mpi_integer4, igq, itag, irbuf9, iend, mpi_integer4, &
	ifq, itag, comm3Dq, istatus, ierr)
	call mpi_sendrecv(rsforce, kpacket9, mpi_real8, igq, itag, rrbuf9, iend, mpi_real8, &
	ifq, itag, comm3Dq, istatus, ierr)

	if(ifq.ne.mpi_proc_null.and.iend.ne.0) then

		jstart = jend + 1
		jend = jend + iend
		jsplits_total = jend

		do j = jstart, jend

			isforce(j) = irbuf9(j - jstart + 1)
			rsforce(j) = rrbuf9(j - jstart + 1)

		enddo

	endif

	!write(*,*) iaq, "s-final", jsplits_total

	return
	end subroutine
