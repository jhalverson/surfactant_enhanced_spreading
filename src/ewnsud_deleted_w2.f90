	subroutine ewnsud_deleted_w2(comm3Dq, iaq, ibq, icq, idq, ieq, ifq, igq)
	
	implicit none
	include "mpif.h"
	include "mybox.cmns"
	
	integer*4 comm3Dq, iaq, ibq, icq, idq, ieq, ifq, igq
	integer*4 itag, istatus(mpi_status_size), ierr
	integer*4 jstart, jend, kstart, kend
	integer*4 j, iend, iend3, iend6
	integer*4 kpacket1, kpacket13, kpacket16
	integer*4 kpacket3, kpacket33, kpacket36
	integer*4 kpacket9, kpacket93, kpacket96
	
	integer*4 irbuf1(kdeletew213), irbuf3(kdeletew233), irbuf9(kdeletew293)
	real*8 rrbuf1(kdeletew216), rrbuf3(kdeletew236), rrbuf9(kdeletew296)
	
	itag = 0
	
	jstart = 1
	jend = jdeletedw2
	jdeletedw2_total = jend
	
	kpacket1 = jdeletedw2
	
	iend = 0
	! send east receive from west
	call mpi_sendrecv(kpacket1, 1, mpi_integer4, ibq, itag, iend, 1, mpi_integer4, &
	icq, itag, comm3Dq, istatus, ierr)
	
	iend3 = 3 * iend
	kpacket13 = 3 * kpacket1
	
	call mpi_sendrecv(iw2ewnsud, kpacket13, mpi_integer4, ibq, itag, irbuf1, iend3, mpi_integer4, &
	icq, itag, comm3Dq, istatus, ierr)
	
	iend6 = 3 * 6 * iend
	kpacket16 = 3 * 6 * kpacket1
	if(iend6.gt.kdeletew216) write(*,*) "ewnsud_deleted_w2 iend6"
	
	call mpi_sendrecv(rw2ewnsud, kpacket16, mpi_real8, ibq, itag, rrbuf1, iend6, mpi_real8, &
	icq, itag, comm3Dq, istatus, ierr)
	
	if(icq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jdeletedw2_total = jend
		
		kstart = 3 * (jstart - 1) + 1
		kend = 3 * jend
		
		do j = kstart, kend
		
			iw2ewnsud(j) = irbuf1(j - kstart + 1)
		
		enddo
		
		kstart = 3 * 6 * (jstart - 1) + 1
		kend = 3 * 6 * jend
		
		do j = kstart, kend
		
			rw2ewnsud(j) = rrbuf1(j - kstart + 1)
		
		enddo
	
	endif
	
	iend = 0
	! send west receive from east
	call mpi_sendrecv(kpacket1, 1, mpi_integer4, icq, itag, iend, 1, mpi_integer4, &
	ibq, itag, comm3Dq, istatus, ierr)
	
	iend3 = 3 * iend
	kpacket13 = 3 * kpacket1
	
	call mpi_sendrecv(iw2ewnsud, kpacket13, mpi_integer4, icq, itag, irbuf1, iend3, mpi_integer4, &
	ibq, itag, comm3Dq, istatus, ierr)
	
	iend6 = 3 * 6 * iend
	kpacket16 = 3 * 6 * kpacket1
	
	call mpi_sendrecv(rw2ewnsud, kpacket16, mpi_real8, icq, itag, rrbuf1, iend6, mpi_real8, &
	ibq, itag, comm3Dq, istatus, ierr)
	
	if(ibq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jdeletedw2_total = jend
		
		kstart = 3 * (jstart - 1) + 1
		kend = 3 * jend
		
		do j = kstart, kend
		
			iw2ewnsud(j) = irbuf1(j - kstart + 1)
			
		enddo
		
		kstart = 3 * 6 * (jstart - 1) + 1
		kend = 3 * 6 * jend
		
		do j = kstart, kend
		
			rw2ewnsud(j) = rrbuf1(j - kstart + 1)
			
		enddo
	
	endif
	
	kpacket3 = jend
	
	iend = 0
	! send north receive from south
	call mpi_sendrecv(kpacket3, 1, mpi_integer4, idq, itag, iend, 1, mpi_integer4, &
	ieq, itag, comm3Dq, istatus, ierr)
	
	iend3 = 3 * iend
	kpacket33 = 3 * kpacket3
	
	call mpi_sendrecv(iw2ewnsud, kpacket33, mpi_integer4, idq, itag, irbuf3, iend3, mpi_integer4, &
	ieq, itag, comm3Dq, istatus, ierr)
	
	iend6 = 3 * 6 * iend
	kpacket36 = 3 * 6 * kpacket3
	if(iend6.gt.kdeletew236) write(*,*) "ewnsud_deleted_w2 iend6"
	
	call mpi_sendrecv(rw2ewnsud, kpacket36, mpi_real8, idq, itag, rrbuf3, iend6, mpi_real8, &
	ieq, itag, comm3Dq, istatus, ierr)
	
	if(ieq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jdeletedw2_total = jend
		
		kstart = 3 * (jstart - 1) + 1
		kend = 3 * jend
		
		do j = kstart, kend
		
			iw2ewnsud(j) = irbuf3(j - kstart + 1)
			
		enddo
		
		kstart = 3 * 6 * (jstart - 1) + 1
		kend = 3 * 6 * jend
		
		do j = kstart, kend
		
			rw2ewnsud(j) = rrbuf3(j - kstart + 1)
			
		enddo
	
	endif
	
	iend = 0
	! send south receive from north
	call mpi_sendrecv(kpacket3, 1, mpi_integer4, ieq, itag, iend, 1, mpi_integer4, &
	idq, itag, comm3Dq, istatus, ierr)
	
	iend3 = 3 * iend
	kpacket33 = 3 * kpacket3
	
	call mpi_sendrecv(iw2ewnsud, kpacket33, mpi_integer4, ieq, itag, irbuf3, iend3, mpi_integer4, &
	idq, itag, comm3Dq, istatus, ierr)
	
	iend6 = 3 * 6 * iend
	kpacket36 = 3 * 6 * kpacket3
	
	call mpi_sendrecv(rw2ewnsud, kpacket36, mpi_real8, ieq, itag, rrbuf3, iend6, mpi_real8, &
	idq, itag, comm3Dq, istatus, ierr)
	
	if(idq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jdeletedw2_total = jend
		
		kstart = 3 * (jstart - 1) + 1
		kend = 3 * jend
		
		do j = kstart, kend
		
			iw2ewnsud(j) = irbuf3(j - kstart + 1)
			
		enddo
		
		kstart = 3 * 6 * (jstart - 1) + 1
		kend = 3 * 6 * jend
		
		do j = kstart, kend
		
			rw2ewnsud(j) = rrbuf3(j - kstart + 1)
			
		enddo
	
	endif
	
	kpacket9 = jend
	
	iend = 0
	! send up receive from down
	call mpi_sendrecv(kpacket9, 1, mpi_integer4, ifq, itag, iend, 1, mpi_integer4, &
	igq, itag, comm3Dq, istatus, ierr)
	
	iend3 = 3 * iend
	kpacket93 = 3 * kpacket9
	
	call mpi_sendrecv(iw2ewnsud, kpacket93, mpi_integer4, ifq, itag, irbuf9, iend3, mpi_integer4, &
	igq, itag, comm3Dq, istatus, ierr)
	
	iend6 = 3 * 6 * iend
	kpacket96 = 3 * 6 * kpacket9
	if(iend6.gt.kdeletew296) write(*,*) "ewnsud_deleted_w2 iend6"
	
	call mpi_sendrecv(rw2ewnsud, kpacket96, mpi_real8, ifq, itag, rrbuf9, iend6, mpi_real8, &
	igq, itag, comm3Dq, istatus, ierr)
	
	if(igq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jdeletedw2_total = jend
		
		kstart = 3 * (jstart - 1) + 1
		kend = 3 * jend
		
		do j = kstart, kend
		
			iw2ewnsud(j) = irbuf9(j - kstart + 1)
			
		enddo
		
		kstart = 3 * 6 * (jstart - 1) + 1
		kend = 3 * 6 * jend
		if(kend.gt.kdeletew236) write(*,*) "ewnsud_deleted_w2 kend"
		
		do j = kstart, kend
		
			rw2ewnsud(j) = rrbuf9(j - kstart + 1)
			
		enddo
	
	endif
	
	iend = 0
	! send up receive from down
	call mpi_sendrecv(kpacket9, 1, mpi_integer4, igq, itag, iend, 1, mpi_integer4, &
	ifq, itag, comm3Dq, istatus, ierr)
	
	iend3 = 3 * iend
	kpacket93 = 3 * kpacket9
	
	call mpi_sendrecv(iw2ewnsud, kpacket93, mpi_integer4, igq, itag, irbuf9, iend3, mpi_integer4, &
	ifq, itag, comm3Dq, istatus, ierr)
	
	iend6 = 3 * 6 * iend
	kpacket96 = 3 * 6 * kpacket9
	
	call mpi_sendrecv(rw2ewnsud, kpacket96, mpi_real8, igq, itag, rrbuf9, iend6, mpi_real8, &
	ifq, itag, comm3Dq, istatus, ierr)
	
	if(ifq.ne.mpi_proc_null.and.iend.ne.0) then
	
		jstart = jend + 1
		jend = jend + iend
		jdeletedw2_total = jend
		
		kstart = 3 * (jstart - 1) + 1
		kend = 3 * jend
		
		do j = kstart, kend
		
			iw2ewnsud(j) = irbuf9(j - kstart + 1)
		
		enddo
		
		kstart = 3 * 6 * (jstart - 1) + 1
		kend = 3 * 6 * jend
		if(kend.gt.kdeletew236) write(*,*) "ewnsud_deleted_w2 kend"
		
		do j = kstart, kend
		
			rw2ewnsud(j) = rrbuf9(j - kstart + 1)
		
		enddo
	
	endif
	
	return
	end subroutine
