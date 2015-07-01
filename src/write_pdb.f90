	subroutine write_pdb(raq, rbq, rcq, iaq, rdq, req, rfq, ibq, icq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq, ibq, icq
	integer*4 i, j, k, iatom, ilimit, jlimit
	integer*4 number_of_files, number_of_full_files, molecules_of_final_file
	real*8 rlimit

	real*8 raq(amxwtr), rbq(amxwtr), rcq(amxwtr)
	real*8 rdq(amxsrf), req(amxsrf), rfq(amxsrf)
	
	character(32) numform
	character(20) flnm
	character(10) substr
	character(10) fldst
	
	numform = '(a6,i5,a5,a4,i6,f12.3,f8.3,f8.3)'
	fldst = 'abcdefghij'
	ilimit = 33333
	jlimit = 99999
	rlimit = 33333.0001D0
	
	write(flnm,'(i7)') icq
	
	if(icq.lt.10) then
	   
	   substr = flnm(7:7)
	   flnm = "000000"//substr
	   flnm = "data"//flnm
	   
	elseif(icq.lt.100) then
	   
	   substr = flnm(6:7)
	   flnm = "00000"//substr
	   flnm = "data"//flnm
	   
	elseif(icq.lt.1000) then
	   
	   substr = flnm(5:7)
	   flnm = "0000"//substr
	   flnm = "data"//flnm
	   
	elseif(icq.lt.10000) then
	   
	   substr = flnm(4:7)
	   flnm = "000"//substr
	   flnm = "data"//flnm
	   
	elseif(icq.lt.100000) then
	   
	   substr = flnm(3:7)
	   flnm = "00"//substr
	   flnm = "data"//flnm
	   
	elseif(icq.lt.1000000) then
	   
	   substr = flnm(2:7)
	   flnm = "0"//substr
	   flnm = "data"//flnm
	   
	else
	   
	   substr = flnm(1:7)
	   flnm = substr
	   flnm = "data"//flnm
	   
	endif
	
	number_of_files = int(iwater / rlimit) + 1
	number_of_full_files = number_of_files - 1
	molecules_of_final_file = iwater - ilimit * number_of_full_files
	
	do k = 1, number_of_full_files
	   
	   flnm = flnm(1:11)//fldst(k:k)
	   flnm = flnm(1:12)//".pdb"
	   
	   open(10, file = flnm)
	   do i = 1, ilimit
	      j = (k - 1) * jlimit + (i - 1) * 3 + 1
	      write(10, numform) "HETATM", (i - 1) * 3 + 1, "O","WAT", i, raq(j + 0)*sigma, rbq(j + 0)*sigma, rcq(j + 0)*sigma
	      write(10, numform) "HETATM", (i - 1) * 3 + 2, "H","WAT", i, raq(j + 1)*sigma, rbq(j + 1)*sigma, rcq(j + 1)*sigma
	      write(10, numform) "HETATM", (i - 1) * 3 + 3, "H","WAT", i, raq(j + 2)*sigma, rbq(j + 2)*sigma, rcq(j + 2)*sigma
	   enddo
	   close(10)
	   
	enddo
	
	flnm = flnm(1:11)//fldst(number_of_files:number_of_files)
	flnm = flnm(1:12)//".pdb"
	
	open(10, file = flnm)
	do i = 1, molecules_of_final_file
	   j = number_of_full_files * jlimit + (i - 1) * 3 + 1
	   write(10, numform) "HETATM", (i - 1) * 3 + 1, "O", "WAT", i, raq(j + 0)*sigma, rbq(j + 0)*sigma, rcq(j + 0)*sigma
	   write(10, numform) "HETATM", (i - 1) * 3 + 2, "H", "WAT", i, raq(j + 1)*sigma, rbq(j + 1)*sigma, rcq(j + 1)*sigma
	   write(10, numform) "HETATM", (i - 1) * 3 + 3, "H", "WAT", i, raq(j + 2)*sigma, rbq(j + 2)*sigma, rcq(j + 2)*sigma
	enddo
	close(10)
	
	flnm = flnm(1:11)//"s.pdb"
	open(10, file = flnm)
	iatom = 0
	do i = 1, ibq
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "H","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "O","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "O","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "O","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "O","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "O","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom,"Si","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "O","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom,"Si","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "O","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom,"Si","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	iatom = iatom + 1
	write(10, numform) "HETATM", iatom, "C","SRF", i, rdq(iatom)*sigma, req(iatom)*sigma, rfq(iatom)*sigma
	enddo
	close(10)

	return
	end subroutine
