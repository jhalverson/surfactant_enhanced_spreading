	subroutine write_force(raq, rbq, rcq, iaq, rdq, req, rfq, ibq, icq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq, ibq, icq
	integer*4 i, iatom

	real*8 raq(amxwtr), rbq(amxwtr), rcq(amxwtr)
	real*8 rdq(amxsrf), req(amxsrf), rfq(amxsrf)
	
	character(32) numform
	character(20) flnm
	character(10) substr
	
	numform = '(a6,i5,a5,a4,i6,f12.3,f8.3,f8.3)'
	
	write(flnm,'(i7)') icq
	
	if(icq.lt.10) then
	   
	   substr = flnm(7:7)
	   flnm = "000000"//substr
	   flnm = "data"//flnm
	   flnm = flnm(1:11)//".pdb"
	   
	elseif(icq.lt.100) then
	   
	   substr = flnm(6:7)
	   flnm = "00000"//substr
	   flnm = "data"//flnm
	   flnm = flnm(1:11)//".pdb"
	   
	elseif(icq.lt.1000) then
	   
	   substr = flnm(5:7)
	   flnm = "0000"//substr
	   flnm = "data"//flnm
	   flnm = flnm(1:11)//".pdb"
	   
	elseif(icq.lt.10000) then
	   
	   substr = flnm(4:7)
	   flnm = "000"//substr
	   flnm = "data"//flnm
	   flnm = flnm(1:11)//".pdb"
	   
	elseif(icq.lt.100000) then
	   
	   substr = flnm(3:7)
	   flnm = "00"//substr
	   flnm = "data"//flnm
	   flnm = flnm(1:11)//".pdb"
	   
	elseif(icq.lt.1000000) then
	   
	   substr = flnm(2:7)
	   flnm = "0"//substr
	   flnm = "data"//flnm
	   flnm = flnm(1:11)//".pdb"
	   
	else
	   
	   substr = flnm(1:7)
	   flnm = substr
	   flnm = "data"//flnm
	   flnm = flnm(1:11)//".pdb"
	   
	endif

	flnm = "f"//flnm
	iatom = 0
	open(10, file = flnm)
        do i = 1, 29 * ibq
        iatom = iatom + 1
        write(10,'(f14.6,f15.6,f15.6)') rdq(iatom), req(iatom), rfq(iatom)
        enddo

	iatom = 0
	do i = 1, 3 * iaq
	iatom = iatom + 1
	write(10,'(f14.6,f15.6,f15.6)') raq(iatom), rbq(iatom), rcq(iatom)
	enddo

	!iatom = 0
	!do i = 1, 29 * ibq
	!iatom = iatom + 1
	!write(10,'(f14.6,f15.6,f15.6)') rdq(iatom), req(iatom), rfq(iatom)
	!enddo
	close(10)

	return
	end subroutine
