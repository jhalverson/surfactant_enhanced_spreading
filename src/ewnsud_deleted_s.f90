	subroutine ewnsud_deleted_s(comm3Dq, iaq, ibq, icq, idq, ieq, ifq, igq)
	
	implicit none
	include "mpif.h"
	include "mybox.cmns"
	
	integer*4 comm3Dq, iaq, ibq, icq, idq, ieq, ifq, igq
	integer*4 itag, istatus(mpi_status_size), ierr
	integer*4 jstart, jend, kstart, kend
	integer*4 j, iend, iend2, iend6
	integer*4 kpacket1, kpacket12, kpacket16
	integer*4 kpacket3, kpacket32, kpacket36
	integer*4 kpacket9, kpacket92, kpacket96
	
	integer*4 irbuf1(kdeletes12), irbuf3(kdeletes32), irbuf9(kdeletes92)
	real*8 rrbuf1(kdeletes16), rrbuf3(kdeletes36), rrbuf9(kdeletes96)
	
	itag = 0
	
	jstart = 1
	jend = jdeleteds
	jdeleteds_total = jend
	
	kpacket1 = jdeleteds
	
	iend = 0
	! send east receive from west
	call mpi_sendrecv(kpacket1, 1, mpi_integer4, ibq, itag, iend, 1, mpi_integer4, &
	icq, itag, comm3Dq, istatus, ierr)
	
	iend2 = 2 * iend
	kpacket12 = 2 * kpacket1
	
	call mpi_sendrecv(isewnsud, kpacket12, mpi_integer4, ibq, itag, irbuf1, iend2, mpi_integer4, &
	icq, itag, comm3Dq, istatus, ierr)
	
	iend6 = 29 * 6 * iend
	kpacket16 = 29 * 6 * kpacket1
	if(iend6.gt.kdeletes16) write(*,*) "ewnsud_deleted_s iend6"
	
	call mpi_sendrecv(rsewnsud, kpacket16, mpi_real8, ibq, itag, rrbuf1, iend6, mpi_real8, &
	icq, itag, comm3Dq, istatus, ierr)
	
	if(icq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jdeleteds_total = jend
		
		kstart = 2 * (jstart - 1) + 1
		kend = 2 * jend
		
		do j = kstart, kend
		
			isewnsud(j) = irbuf1(j - kstart + 1)
		
		enddo
		
		kstart = 29 * 6 * (jstart - 1) + 1
		kend = 29 * 6 * jend
		
		do j = kstart, kend
		
			rsewnsud(j) = rrbuf1(j - kstart + 1)
		
		enddo
	
	endif
	
	iend = 0
	! send west receive from east
	call mpi_sendrecv(kpacket1, 1, mpi_integer4, icq, itag, iend, 1, mpi_integer4, &
	ibq, itag, comm3Dq, istatus, ierr)
	
	iend2 = 2 * iend
	kpacket12 = 2 * kpacket1
	
	call mpi_sendrecv(isewnsud, kpacket12, mpi_integer4, icq, itag, irbuf1, iend2, mpi_integer4, &
	ibq, itag, comm3Dq, istatus, ierr)
	
	iend6 = 29 * 6 * iend
	kpacket16 = 29 * 6 * kpacket1
	
	call mpi_sendrecv(rsewnsud, kpacket16, mpi_real8, icq, itag, rrbuf1, iend6, mpi_real8, &
	ibq, itag, comm3Dq, istatus, ierr)
	
	if(ibq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jdeleteds_total = jend
		
		kstart = 2 * (jstart - 1) + 1
		kend = 2 * jend
		
		do j = kstart, kend
		
			isewnsud(j) = irbuf1(j - kstart + 1)
			
		enddo
		
		kstart = 29 * 6 * (jstart - 1) + 1
		kend = 29 * 6 * jend
		
		do j = kstart, kend
		
			rsewnsud(j) = rrbuf1(j - kstart + 1)
			
		enddo
	
	endif
	
	kpacket3 = jend
	
	iend = 0
	! send north receive from south
	call mpi_sendrecv(kpacket3, 1, mpi_integer4, idq, itag, iend, 1, mpi_integer4, &
	ieq, itag, comm3Dq, istatus, ierr)
	
	iend2 = 2 * iend
	kpacket32 = 2 * kpacket3
	
	call mpi_sendrecv(isewnsud, kpacket32, mpi_integer4, idq, itag, irbuf3, iend2, mpi_integer4, &
	ieq, itag, comm3Dq, istatus, ierr)
	
	iend6 = 29 * 6 * iend
	kpacket36 = 29 * 6 * kpacket3
	if(iend6.gt.kdeletes36) write(*,*) "ewnsud_deleted_s iend6"
	
	call mpi_sendrecv(rsewnsud, kpacket36, mpi_real8, idq, itag, rrbuf3, iend6, mpi_real8, &
	ieq, itag, comm3Dq, istatus, ierr)
	
	if(ieq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jdeleteds_total = jend
		
		kstart = 2 * (jstart - 1) + 1
		kend = 2 * jend
		
		do j = kstart, kend
		
			isewnsud(j) = irbuf3(j - kstart + 1)
			
		enddo
		
		kstart = 29 * 6 * (jstart - 1) + 1
		kend = 29 * 6 * jend
		
		do j = kstart, kend
		
			rsewnsud(j) = rrbuf3(j - kstart + 1)
			
		enddo
	
	endif
	
	iend = 0
	! send south receive from north
	call mpi_sendrecv(kpacket3, 1, mpi_integer4, ieq, itag, iend, 1, mpi_integer4, &
	idq, itag, comm3Dq, istatus, ierr)
	
	iend2 = 2 * iend
	kpacket32 = 2 * kpacket3
	
	call mpi_sendrecv(isewnsud, kpacket32, mpi_integer4, ieq, itag, irbuf3, iend2, mpi_integer4, &
	idq, itag, comm3Dq, istatus, ierr)
	
	iend6 = 29 * 6 * iend
	kpacket36 = 29 * 6 * kpacket3
	
	call mpi_sendrecv(rsewnsud, kpacket36, mpi_real8, ieq, itag, rrbuf3, iend6, mpi_real8, &
	idq, itag, comm3Dq, istatus, ierr)
	
	if(idq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jdeleteds_total = jend
		
		kstart = 2 * (jstart - 1) + 1
		kend = 2 * jend
		
		do j = kstart, kend
		
			isewnsud(j) = irbuf3(j - kstart + 1)
			
		enddo
		
		kstart = 29 * 6 * (jstart - 1) + 1
		kend = 29 * 6 * jend
		
		do j = kstart, kend
		
			rsewnsud(j) = rrbuf3(j - kstart + 1)
			
		enddo
	
	endif
	
	kpacket9 = jend
	
	iend = 0
	! send up receive from down
	call mpi_sendrecv(kpacket9, 1, mpi_integer4, ifq, itag, iend, 1, mpi_integer4, &
	igq, itag, comm3Dq, istatus, ierr)
	
	iend2 = 2 * iend
	kpacket92 = 2 * kpacket9
	
	call mpi_sendrecv(isewnsud, kpacket92, mpi_integer4, ifq, itag, irbuf9, iend2, mpi_integer4, &
	igq, itag, comm3Dq, istatus, ierr)
	
	iend6 = 29 * 6 * iend
	kpacket96 = 29 * 6 * kpacket9
	if(iend6.gt.kdeletes96) write(*,*) "ewnsud_deleted_s iend6"
	
	call mpi_sendrecv(rsewnsud, kpacket96, mpi_real8, ifq, itag, rrbuf9, iend6, mpi_real8, &
	igq, itag, comm3Dq, istatus, ierr)
	
	if(igq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jdeleteds_total = jend
		
		kstart = 2 * (jstart - 1) + 1
		kend = 2 * jend
		
		do j = kstart, kend
		
			isewnsud(j) = irbuf9(j - kstart + 1)
			
		enddo
		
		kstart = 29 * 6 * (jstart - 1) + 1
		kend = 29 * 6 * jend
		if(kend.gt.kdeletes296) write(*,*) "ewnsud_deleted_s kend"
		
		do j = kstart, kend
		
			rsewnsud(j) = rrbuf9(j - kstart + 1)
			
		enddo
	
	endif
	
	iend = 0
	! send up receive from down
	call mpi_sendrecv(kpacket9, 1, mpi_integer4, igq, itag, iend, 1, mpi_integer4, &
	ifq, itag, comm3Dq, istatus, ierr)
	
	iend2 = 2 * iend
	kpacket92 = 2 * kpacket9
	
	call mpi_sendrecv(isewnsud, kpacket92, mpi_integer4, igq, itag, irbuf9, iend2, mpi_integer4, &
	ifq, itag, comm3Dq, istatus, ierr)
	
	iend6 = 29 * 6 * iend
	kpacket96 = 29 * 6 * kpacket9
	
	call mpi_sendrecv(rsewnsud, kpacket96, mpi_real8, igq, itag, rrbuf9, iend6, mpi_real8, &
	ifq, itag, comm3Dq, istatus, ierr)
	
	if(ifq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jdeleteds_total = jend
		
		kstart = 2 * (jstart - 1) + 1
		kend = 2 * jend
		
		do j = kstart, kend
		
			isewnsud(j) = irbuf9(j - kstart + 1)
		
		enddo
		
		kstart = 29 * 6 * (jstart - 1) + 1
		kend = 29 * 6 * jend
		if(kend.gt.kdeletes296) write(*,*) "ewnsud_deleted_s kend"
		
		do j = kstart, kend
		
			rsewnsud(j) = rrbuf9(j - kstart + 1)
		
		enddo
	
	endif
	
	return
	end subroutine
