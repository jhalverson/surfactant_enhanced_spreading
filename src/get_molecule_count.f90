	subroutine get_molecule_count(iaq)

	implicit none
	include "mybox.cmns"
	include "mpif.h"

	integer*4 iaq
	integer*4 j, icnt, icntotal, ierr
	integer*4 ib1000, ib64to1000
	integer*4 ib512, ib64to512
	integer*4 ib216, ib64to216

	integer*4 ictabcw1, ictabcw2, ictabcs, itotatms
	integer*4 ictabcw1total, ictabcw2total, ictabcstotal

	icnt = 0
	itotatms = 0
	
	ictabcw1=0
	ictabcw2=0
	ictabcs=0

	do j = 1, 64

		ib216 = ib64to216(j)
		icnt = icnt + ictw1(ib216)
		ictabcw1 = ictabcw1 + ictw1(ib216)
		itotatms = itotatms + 3 * ictw1(ib216)
		!write(*,*) iaq, j, ictw1(ib216), "w1"

	enddo

	do j = 1, 64

		ib512 = ib64to512(j)
		icnt = icnt + ictw2(ib512)
		ictabcw2 = ictabcw2 + ictw2(ib512)
		itotatms = itotatms + 3 * ictw2(ib512)
		!if(iaq.eq.13) write(*,*) iaq, j, ictw2(ib216), "w2"

	enddo

	do j = 1, 64

		ib1000 = ib64to1000(j)
		icnt = icnt + icts(ib1000)
		ictabcs = ictabcs + icts(ib1000)
		itotatms = itotatms + 29 * icts(ib1000)
		!if(iaq.eq.13) write(*,*) iaq, j, icts(ib1000), "s"

	enddo

	!write(*,'(i5,i5,i5,i5)') iaq, ictabcw1, ictabcw2, ictabcs
	!write(*,*) iaq, itotatms

	call mpi_reduce(icnt, icntotal, 1, mpi_integer4, mpi_sum, 0, mpi_comm_world, ierr)
	if(iaq.eq.0.and.icntotal.ne.iwater+isurfactant) write(*,*) "MOLECULES LOST:", icntotal, "OF", iwater+isurfactant 
	call mpi_reduce(ictabcw1, ictabcw1total, 1, mpi_integer4, mpi_sum, 0, mpi_comm_world, ierr)
	call mpi_reduce(ictabcw2, ictabcw2total, 1, mpi_integer4, mpi_sum, 0, mpi_comm_world, ierr)
	call mpi_reduce(ictabcs, ictabcstotal, 1, mpi_integer4, mpi_sum, 0, mpi_comm_world, ierr)
	!if(iaq.eq.0) write(*,*) "molecules", ictabcw1total, ictabcw2total, ictabcstotal
	return
	end subroutine
