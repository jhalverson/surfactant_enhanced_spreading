	subroutine force_to_root(iaq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"

	integer*4 iaq
	integer*4 ict, i, j, k, m, ierr
	integer*4 ib1000, ib64to1000
	integer*4 ib512, ib64to512
	integer*4 ib216, ib64to216

	real*8 rcnfgx(amx4), rcnfgy(amx4), rcnfgz(amx4)

	ict = 0

	do j = 1, 64

		ib216 = ib64to216(j)

		do i = 1, ictw1(ib216)

			do k = 1, 3

				ict = ict + 1
				m = (i - 1) * 3 + k
				rcnfgx(ict) = fxw1(m, ib216)
				rcnfgy(ict) = fyw1(m, ib216)
				rcnfgz(ict) = fzw1(m, ib216)

			enddo

		enddo

	enddo

	do j = 1, 64

		ib512 = ib64to512(j)

		do i = 1, ictw2(ib512)

			do k = 1, 3

				ict = ict + 1
				m = (i - 1) * 3 + k
				rcnfgx(ict) = fxw2(m, ib512)
				rcnfgy(ict) = fyw2(m, ib512)
				rcnfgz(ict) = fzw2(m, ib512)

			enddo

		enddo

	enddo

	call mpi_send(ict, 1, mpi_integer4, 0, iaq, mpi_comm_world, ierr)
	call mpi_send(rcnfgx, ict, mpi_real8, 0, iaq, mpi_comm_world, ierr)
	call mpi_send(rcnfgy, ict, mpi_real8, 0, iaq, mpi_comm_world, ierr)
	call mpi_send(rcnfgz, ict, mpi_real8, 0, iaq, mpi_comm_world, ierr)

	ict = 0

	do j = 1, 64

		ib1000 = ib64to1000(j)

		do i = 1, icts(ib1000)

			do k = 1, 29

				ict = ict + 1
				m = (i - 1) * 29 + k
				rcnfgx(ict) = fxs(m, ib1000)
				rcnfgy(ict) = fys(m, ib1000)
				rcnfgz(ict) = fzs(m, ib1000)

			enddo

		enddo

	enddo

	call mpi_send(ict, 1, mpi_integer4, 0, iaq, mpi_comm_world, ierr)
	call mpi_send(rcnfgx, ict, mpi_real8, 0, iaq, mpi_comm_world, ierr)
	call mpi_send(rcnfgy, ict, mpi_real8, 0, iaq, mpi_comm_world, ierr)
	call mpi_send(rcnfgz, ict, mpi_real8, 0, iaq, mpi_comm_world, ierr)

	return
	end subroutine
