	subroutine ewnsud_Onm2(comm3Dq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"

	integer*4 comm3Dq, ierr

	call mpi_allgather(cO2, lp, mpi_complex16, cO2all, lp, mpi_complex16, comm3Dq, ierr)

	return
	end subroutine
