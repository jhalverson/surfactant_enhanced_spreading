	subroutine force_from_root(iaq, ibq, icq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"

	integer*4 iaq, ibq, icq
	integer*4 j, i, k, m
	integer*4 iatomw, iatoms, ictww, ictss, iwtrmle, isrfmle
	integer*4 ib1000, ib64to1000, ib512, ib64to512, ib216, ib64to216
	integer*4 istatus(mpi_status_size), ierr

	real*8 xt(amxwtr), yt(amxwtr), zt(amxwtr)
	real*8 xsurf(amxsrf), ysurf(amxsrf), zsurf(amxsrf)
	real*8 rcgx(amx4), rcgy(amx4), rcgz(amx4)

	iatomw = 0
	iatoms = 0

	iwtrmle = 0
	isrfmle = 0

	do j = 1, 64

		ib216 = ib64to216(j)

		do i = 1, ictw1(ib216)

			iwtrmle = iwtrmle + 1

			do k = 1, 3

				m = (i - 1) * 3 + k
				iatomw = iatomw + 1
				xt(iatomw) = fxw1(m, ib216)
				yt(iatomw) = fyw1(m, ib216)
				zt(iatomw) = fzw1(m, ib216)

			enddo

		enddo

	enddo

	do j = 1, 64

		ib512 = ib64to512(j)

		do i = 1, ictw2(ib512)

			iwtrmle = iwtrmle + 1

			do k = 1, 3

				m = (i - 1) * 3 + k
				iatomw = iatomw + 1
				xt(iatomw) = fxw2(m, ib512)
				yt(iatomw) = fyw2(m, ib512)
				zt(iatomw) = fzw2(m, ib512)

			enddo

		enddo

	enddo

	do j = 1, 64

		ib1000 = ib64to1000(j)

		do i = 1, icts(ib1000)

			isrfmle = isrfmle + 1

			do k = 1, 29

				m = (i - 1) * 29 + k
				iatoms = iatoms + 1
				xsurf(iatoms) = fxs(m, ib1000)
				ysurf(iatoms) = fys(m, ib1000)
				zsurf(iatoms) = fzs(m, ib1000)

			enddo

		enddo

	enddo

	do i = 1, ibq - 1

		call mpi_recv(ictww, 1, mpi_integer4, i, i, mpi_comm_world, istatus, ierr)
		call mpi_recv(rcgx, ictww, mpi_real8, i, i, mpi_comm_world, istatus, ierr)
		call mpi_recv(rcgy, ictww, mpi_real8, i, i, mpi_comm_world, istatus, ierr)
		call mpi_recv(rcgz, ictww, mpi_real8, i, i, mpi_comm_world, istatus, ierr)

		do j = 1, ictww

			iatomw = iatomw + 1
			xt(iatomw) = rcgx(j)
			yt(iatomw) = rcgy(j)
			zt(iatomw) = rcgz(j)

		enddo

		iwtrmle = iwtrmle + ictww / 3

		call mpi_recv(ictss, 1, mpi_integer4, i, i, mpi_comm_world, istatus, ierr)
		call mpi_recv(rcgx, ictss, mpi_real8, i, i, mpi_comm_world, istatus, ierr)
		call mpi_recv(rcgy, ictss, mpi_real8, i, i, mpi_comm_world, istatus, ierr)
		call mpi_recv(rcgz, ictss, mpi_real8, i, i, mpi_comm_world, istatus, ierr)

		do j = 1, ictss

			iatoms = iatoms + 1
			xsurf(iatoms) = rcgx(j)
			ysurf(iatoms) = rcgy(j)
			zsurf(iatoms) = rcgz(j)

		enddo

		isrfmle = isrfmle + ictss / 29

	enddo
	!write(*,*) "count", iwtrmle, isrfmle
	call write_force(xt, yt, zt, iwtrmle, xsurf, ysurf, zsurf, isrfmle, icq)

	return
	end subroutine
