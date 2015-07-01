	subroutine init_box_dimensions()

	implicit none
	include "mybox.cmns"

	ibdx3 = 2 * ibdx2
	ibdy3 = 2 * ibdy2
	ibdz3 = 2 * ibdz2

	ibdx4 = 2 * ibdx3
	ibdy4 = 2 * ibdy3
	ibdz4 = 2 * ibdz3

	rsdx = ibdx4 * rs
	rsdy = ibdy4 * rs
	rsdz = ibdz4 * rs

	rhsx = rsdx / 2.0D0
	rhsy = rsdy / 2.0D0
	rhsz = rsdz / 2.0D0

	return
	end subroutine


	subroutine get_box_center(comm3Dq, iaq, ibq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"

	integer*4 comm3Dq, iaq, ibq(3)
	integer*4 i, j, k, l, ierr
	real*8 ribdx2, ribdy2, ribdz2

	ribdx2 = ibdx2
	ribdy2 = ibdy2
	ribdz2 = ibdz2

	rcen2xtmp = rsdx * (ibq(3) + 1) / ribdx2 - rhsx - rsdx / ribdx2 / 2.0D0
	rcen2ytmp = rsdy * (ibq(2) + 1) / ribdy2 - rhsy - rsdy / ribdy2 / 2.0D0
	rcen2ztmp = rsdz * (ibq(1) + 1) / ribdz2 - rhsz - rsdz / ribdz2 / 2.0D0
	!write(*,'(i4, f8.4,f9.4,f9.4)') iaq + 1, rcen2xtmp, rcen2ytmp, rcen2ztmp
	call mpi_allgather(rcen2xtmp, 1, mpi_real8, rcen2x, 1, mpi_real8, comm3Dq, ierr)
	call mpi_allgather(rcen2ytmp, 1, mpi_real8, rcen2y, 1, mpi_real8, comm3Dq, ierr)
	call mpi_allgather(rcen2ztmp, 1, mpi_real8, rcen2z, 1, mpi_real8, comm3Dq, ierr)

	do k = 1, ibdz2

		do j = 1, ibdy2

			do i = 1, ibdx2

				l = i + (j - 1) * ibdx2 + (k - 1) * ibdx2 * ibdy2

			enddo

		enddo

	enddo

	do k = 1, 2

		do j = 1, 2

			do i = 1, 2

				l = i + (j - 1) * 2 + (k - 1) * 2 * 2
				rcen3x(l) = 4.0D0 * rs * i / 2.0D0 - 2.0D0 * rs - 4.0D0 * rs / 4.0D0 + rcen2xtmp
				rcen3y(l) = 4.0D0 * rs * j / 2.0D0 - 2.0D0 * rs - 4.0D0 * rs / 4.0D0 + rcen2ytmp
				rcen3z(l) = 4.0D0 * rs * k / 2.0D0 - 2.0D0 * rs - 4.0D0 * rs / 4.0D0 + rcen2ztmp

			enddo

		enddo

	enddo

	do k = 1, 4

		do j = 1, 4

			do i = 1, 4

				l = i + (j - 1) * 4 + (k - 1) * 4 * 4
				rcen4x(l) = 4.0D0 * rs * i / 4.0D0 - 2.0D0 * rs - 4.0D0 * rs / 8.0D0 + rcen2xtmp
				rcen4y(l) = 4.0D0 * rs * j / 4.0D0 - 2.0D0 * rs - 4.0D0 * rs / 8.0D0 + rcen2ytmp
				rcen4z(l) = 4.0D0 * rs * k / 4.0D0 - 2.0D0 * rs - 4.0D0 * rs / 8.0D0 + rcen2ztmp
				!if(iaq.eq.0) write(*,'(i6,f8.3,f8.3,f8.3)') l,rcen4x(l),rcen4y(l),rcen4z(l)
			enddo

		enddo

	enddo

	! get centers for the 216 L3 boxes of Onm3

	do k = 1, 6

		do j = 1, 6

			do i = 1, 6

				l = i + (j - 1) * 6 + (k - 1) * 6 * 6
				rcen3xLjk3(l) = 12.0D0 * rs * i / 6.0D0 - 6.0D0 * rs - 12.0D0 * rs / 12.0D0 + rcen2xtmp
				rcen3yLjk3(l) = 12.0D0 * rs * j / 6.0D0 - 6.0D0 * rs - 12.0D0 * rs / 12.0D0 + rcen2ytmp
				rcen3zLjk3(l) = 12.0D0 * rs * k / 6.0D0 - 6.0D0 * rs - 12.0D0 * rs / 12.0D0 + rcen2ztmp

			enddo

		enddo

	enddo

	! get centers for the 512 L4 boxes of Onm4

	do k = 1, 8

		do j = 1, 8

			do i = 1, 8

				l = i + (j - 1) * 8 + (k - 1) * 8 * 8
				rcen4xLjk4(l) = 8.0D0 * rs * i / 8.0D0 - 8.0D0 * rs / 2.0D0 - 8.0D0 * rs / 16.0D0 + rcen2xtmp
				rcen4yLjk4(l) = 8.0D0 * rs * j / 8.0D0 - 8.0D0 * rs / 2.0D0 - 8.0D0 * rs / 16.0D0 + rcen2ytmp
				rcen4zLjk4(l) = 8.0D0 * rs * k / 8.0D0 - 8.0D0 * rs / 2.0D0 - 8.0D0 * rs / 16.0D0 + rcen2ztmp
				!if(iaq.eq.0) write(*,'(i6,f8.3,f8.3,f8.3)') l,rcen4xLjk4(l),rcen4yLjk4(l),rcen4zLjk4(l)
			enddo

		enddo

	enddo
	
	! get centers for the 512 w2 boxes
	
	do k = 1, 8

		do j = 1, 8

			do i = 1, 8

				l = i + (j - 1) * 8 + (k - 1) * 8 * 8
				rcen512x(l) = 8.0D0 * rs * i / 8.0D0 - 8.0D0 * rs / 2.0D0 - 8.0D0 * rs / 16.0D0 + rcen2xtmp
				rcen512y(l) = 8.0D0 * rs * j / 8.0D0 - 8.0D0 * rs / 2.0D0 - 8.0D0 * rs / 16.0D0 + rcen2ytmp
				rcen512z(l) = 8.0D0 * rs * k / 8.0D0 - 8.0D0 * rs / 2.0D0 - 8.0D0 * rs / 16.0D0 + rcen2ztmp

			enddo

		enddo

	enddo
	
	! get centers for the 1000 surfactant boxes
	
	do k = 1, 10

		do j = 1, 10

			do i = 1, 10

				l = i + (j - 1) * 10 + (k - 1) * 10 * 10
				rcen1000x(l) = 10.0D0 * rs * i / 10.0D0 - 10.0D0 * rs / 2.0D0 - 10.0D0 * rs / 20.0D0 + rcen2xtmp
				rcen1000y(l) = 10.0D0 * rs * j / 10.0D0 - 10.0D0 * rs / 2.0D0 - 10.0D0 * rs / 20.0D0 + rcen2ytmp
				rcen1000z(l) = 10.0D0 * rs * k / 10.0D0 - 10.0D0 * rs / 2.0D0 - 10.0D0 * rs / 20.0D0 + rcen2ztmp

			enddo

		enddo

	enddo

	return
	end subroutine
