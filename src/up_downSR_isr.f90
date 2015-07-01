	subroutine up_downSR_isr(comm3Dq, iaq, ibq, icq, idq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"
	
	integer*4 comm3Dq, iaq, ibq, icq, idq
	integer*4 j, k, l, m, n
	integer*4 kd, ld
	integer*4 ib1000, ib512, ib216, ictmolecule
	integer*4 i3l, ictcrd
	integer*4 i3ld, ictcrdd
	integer*4 ireq(8), istatus8(mpi_status_size, 8)
	integer*4 itag, ierr, istatus(mpi_status_size)

	integer*4 ictbuf(464), ictrbuf(464)
	integer*4 ictbufd(464), ictrbufd(464)
	real*8 rcrdbuf(kudRws), rcrdrbuf(kudRws)
	real*8 rcrdbufd(kudRws), rcrdrbufd(kudRws)

	! load arrays for sending

	k = 0
	l = 0
	do ib1000 = 401, 500

		k = k + 1
		ictbuf(k) = icts(ib1000)
		ictmolecule = icts(ib1000)

		do j = 1, ictmolecule

			do m = 1, 29

				n = (j - 1) * 29 + m
				rcrdbuf(l + 1) = rxs(n, ib1000)
				rcrdbuf(l + 2) = rys(n, ib1000)
				rcrdbuf(l + 3) = rzs(n, ib1000)
				l = l + 3

			enddo

		enddo

	enddo

	do ib1000 = 501, 600

		k = k + 1
		ictbuf(k) = icts(ib1000)
		ictmolecule = icts(ib1000)

		do j = 1, ictmolecule

			do m = 1, 29

				n = (j - 1) * 29 + m
				rcrdbuf(l + 1) = rxs(n, ib1000)
				rcrdbuf(l + 2) = rys(n, ib1000)
				rcrdbuf(l + 3) = rzs(n, ib1000)
				l = l + 3

			enddo

		enddo

	enddo

	do ib512 = 257, 320

		k = k + 1
		ictbuf(k) = ictw2(ib512)
		ictmolecule = ictw2(ib512)

		do j = 1, ictmolecule

			do m = 1, 3

				n = (j - 1) * 3 + m
				rcrdbuf(l + 1) = rxw2(n, ib512)
				rcrdbuf(l + 2) = ryw2(n, ib512)
				rcrdbuf(l + 3) = rzw2(n, ib512)
				l = l + 3

			enddo

		enddo

	enddo

	do ib1000 = 601, 700

		k = k + 1
		ictbuf(k) = icts(ib1000)
		ictmolecule = icts(ib1000)

		do j = 1, ictmolecule

			do m = 1, 29

				n = (j - 1) * 29 + m
				rcrdbuf(l + 1) = rxs(n, ib1000)
				rcrdbuf(l + 2) = rys(n, ib1000)
				rcrdbuf(l + 3) = rzs(n, ib1000)
				l = l + 3

			enddo

		enddo

	enddo

	do ib512 = 321, 384

		k = k + 1
		ictbuf(k) = ictw2(ib512)
		ictmolecule = ictw2(ib512)

		do j = 1, ictmolecule

			do m = 1, 3

				n = (j - 1) * 3 + m
				rcrdbuf(l + 1) = rxw2(n, ib512)
				rcrdbuf(l + 2) = ryw2(n, ib512)
				rcrdbuf(l + 3) = rzw2(n, ib512)
				l = l + 3

			enddo

		enddo

	enddo

	do ib216 = 145, 180

		k = k + 1
		ictbuf(k) = ictw1(ib216)
		ictmolecule = ictw1(ib216)

		do j = 1, ictmolecule

			do m = 1, 3

				n = (j - 1) * 3 + m
				rcrdbuf(l + 1) = rxw1(n, ib216)
				rcrdbuf(l + 2) = ryw1(n, ib216)
				rcrdbuf(l + 3) = rzw1(n, ib216)
				l = l + 3

			enddo

		enddo

	enddo

	! check
	if(l.gt.kudRws) write(*,*) "upSR l kudRws"
	if(k.ne.464) write(*,*) "upSR k load"

	i3l = l
	
	! load arrays for sending

	kd = 0
	ld = 0
	do ib1000 = 501, 600

		kd = kd + 1
		ictbufd(kd) = icts(ib1000)
		ictmolecule = icts(ib1000)

		do j = 1, ictmolecule

			do m = 1, 29

				n = (j - 1) * 29 + m
				rcrdbufd(ld + 1) = rxs(n, ib1000)
				rcrdbufd(ld + 2) = rys(n, ib1000)
				rcrdbufd(ld + 3) = rzs(n, ib1000)
				ld = ld + 3

			enddo

		enddo

	enddo

	do ib1000 = 401, 500

		kd = kd + 1
		ictbufd(kd) = icts(ib1000)
		ictmolecule = icts(ib1000)

		do j = 1, ictmolecule

			do m = 1, 29

				n = (j - 1) * 29 + m
				rcrdbufd(ld + 1) = rxs(n, ib1000)
				rcrdbufd(ld + 2) = rys(n, ib1000)
				rcrdbufd(ld + 3) = rzs(n, ib1000)
				ld = ld + 3

			enddo

		enddo

	enddo

	do ib512 = 193, 256

		kd = kd + 1
		ictbufd(kd) = ictw2(ib512)
		ictmolecule = ictw2(ib512)

		do j = 1, ictmolecule

			do m = 1, 3

				n = (j - 1) * 3 + m
				rcrdbufd(ld + 1) = rxw2(n, ib512)
				rcrdbufd(ld + 2) = ryw2(n, ib512)
				rcrdbufd(ld + 3) = rzw2(n, ib512)
				ld = ld + 3

			enddo

		enddo

	enddo

	do ib1000 = 301, 400

		kd = kd + 1
		ictbufd(kd) = icts(ib1000)
		ictmolecule = icts(ib1000)

		do j = 1, ictmolecule

			do m = 1, 29

				n = (j - 1) * 29 + m
				rcrdbufd(ld + 1) = rxs(n, ib1000)
				rcrdbufd(ld + 2) = rys(n, ib1000)
				rcrdbufd(ld + 3) = rzs(n, ib1000)
				ld = ld + 3

			enddo

		enddo

	enddo

	do ib512 = 129, 192

		kd = kd + 1
		ictbufd(kd) = ictw2(ib512)
		ictmolecule = ictw2(ib512)

		do j = 1, ictmolecule

			do m = 1, 3

				n = (j - 1) * 3 + m
				rcrdbufd(ld + 1) = rxw2(n, ib512)
				rcrdbufd(ld + 2) = ryw2(n, ib512)
				rcrdbufd(ld + 3) = rzw2(n, ib512)
				ld = ld + 3

			enddo

		enddo

	enddo

	do ib216 = 37, 72

		kd = kd + 1
		ictbufd(kd) = ictw1(ib216)
		ictmolecule = ictw1(ib216)

		do j = 1, ictmolecule

			do m = 1, 3

				n = (j - 1) * 3 + m
				rcrdbufd(ld + 1) = rxw1(n, ib216)
				rcrdbufd(ld + 2) = ryw1(n, ib216)
				rcrdbufd(ld + 3) = rzw1(n, ib216)
				ld = ld + 3

			enddo

		enddo

	enddo

	! check
	if(ld.gt.kudRws) write(*,*) "downSR ld kudRws"
	if(kd.ne.464) write(*,*) "downSR kd load"
	
	i3ld = ld

	itag = 0

	! send up receive from down

	call mpi_sendrecv(i3l, 1, mpi_integer4, ibq, itag, ictcrd, 1, mpi_integer4, &
	icq, itag, comm3Dq, istatus, ierr)
	call mpi_sendrecv(i3ld, 1, mpi_integer4, icq, itag, ictcrdd, 1, mpi_integer4, &
	ibq, itag, comm3Dq, istatus, ierr)
	
	call mpi_isend(ictbuf, 464, mpi_integer4, ibq, itag, comm3Dq, ireq(1), ierr)
	call mpi_isend(rcrdbuf, i3l, mpi_real8, ibq, itag, comm3Dq, ireq(2), ierr)
	call mpi_isend(ictbufd, 464, mpi_integer4, icq, itag, comm3Dq, ireq(3), ierr)
	call mpi_isend(rcrdbufd, i3ld, mpi_real8, icq, itag, comm3Dq, ireq(4), ierr)
	
	call mpi_irecv(ictrbuf, 464, mpi_integer4, icq, itag, comm3Dq, ireq(5), ierr)
	call mpi_irecv(rcrdrbuf, ictcrd, mpi_real8, icq, itag, comm3Dq, ireq(6), ierr)
	call mpi_irecv(ictrbufd, 464, mpi_integer4, ibq, itag, comm3Dq, ireq(7), ierr)
	call mpi_irecv(rcrdrbufd, ictcrdd, mpi_real8, ibq, itag, comm3Dq, ireq(8), ierr)
	
	call evalw1w1(iaq)
	call evalw1w2(iaq)
	call evalw1s(iaq)
	call evalw2w2(iaq)
	call evalw2s(iaq)
	call evalss(iaq)
	
	call eval_valence(iaq)
	call eval_dihedral(iaq)
	call eval_dihedral_berendsen(iaq)
	call eval_surfactant_internal(iaq)
	call eval_subtract_fma(iaq)
	
	call eval_wall()
	call eval_graphite(iaq, idq)
	
	call mpi_waitall(8, ireq, istatus8, ierr)

	! unload arrays if received

	if(icq.ne.mpi_proc_null) then

		k = 0
		l = 0
		do ib1000 = 1, 100

			k = k + 1
			icts(ib1000) = ictrbuf(k)
			ictmolecule = icts(ib1000)

			do j = 1, ictmolecule

				do m = 1, 29

					n = (j - 1) * 29 + m
					rxs(n, ib1000) = rcrdrbuf(l + 1)
					rys(n, ib1000) = rcrdrbuf(l + 2)
					rzs(n, ib1000) = rcrdrbuf(l + 3)
					l = l + 3

				enddo

			enddo

		enddo

		do ib1000 = 101, 200

			k = k + 1
			icts(ib1000) = ictrbuf(k)
			ictmolecule = icts(ib1000)

			do j = 1, ictmolecule

				do m = 1, 29

					n = (j - 1) * 29 + m
					rxs(n, ib1000) = rcrdrbuf(l + 1)
					rys(n, ib1000) = rcrdrbuf(l + 2)
					rzs(n, ib1000) = rcrdrbuf(l + 3)
					l = l + 3

				enddo

			enddo

		enddo

		do ib512 = 1, 64

			k = k + 1
			ictw2(ib512) = ictrbuf(k)
			ictmolecule = ictw2(ib512)

			do j = 1, ictmolecule

				do m = 1, 3

					n = (j - 1) * 3 + m
					rxw2(n, ib512) = rcrdrbuf(l + 1)
					ryw2(n, ib512) = rcrdrbuf(l + 2)
					rzw2(n, ib512) = rcrdrbuf(l + 3)
					l = l + 3

				enddo

			enddo

		enddo

		do ib1000 = 201, 300 

			k = k + 1
			icts(ib1000) = ictrbuf(k)
			ictmolecule = icts(ib1000)

			do j = 1, ictmolecule

				do m = 1, 29

					n = (j - 1) * 29 + m
					rxs(n, ib1000) = rcrdrbuf(l + 1)
					rys(n, ib1000) = rcrdrbuf(l + 2)
					rzs(n, ib1000) = rcrdrbuf(l + 3)
					l = l + 3

				enddo

			enddo

		enddo

		do ib512 = 65, 128

			k = k + 1
			ictw2(ib512) = ictrbuf(k)
			ictmolecule = ictw2(ib512)

			do j = 1, ictmolecule

				do m = 1, 3

					n = (j - 1) * 3 + m
					rxw2(n, ib512) = rcrdrbuf(l + 1)
					ryw2(n, ib512) = rcrdrbuf(l + 2)
					rzw2(n, ib512) = rcrdrbuf(l + 3)
					l = l + 3

				enddo

			enddo

		enddo

		do ib216 = 1, 36

			k = k + 1
			ictw1(ib216) = ictrbuf(k)
			ictmolecule = ictw1(ib216)

			do j = 1, ictmolecule

				do m = 1, 3

					n = (j - 1) * 3 + m
					rxw1(n, ib216) = rcrdrbuf(l + 1)
					ryw1(n, ib216) = rcrdrbuf(l + 2)
					rzw1(n, ib216) = rcrdrbuf(l + 3)
					l = l + 3

				enddo

			enddo

		enddo

		if(l.ne.ictcrd) write(*,*) "upSR l"
		if(k.ne.464) write(*,*) "upSR k unload"

	endif
	
	! unload arrays if received

	if(ibq.ne.mpi_proc_null) then

		kd = 0
		ld = 0
		do ib1000 = 901, 1000

			kd = kd + 1
			icts(ib1000) = ictrbufd(kd)
			ictmolecule = icts(ib1000)

			do j = 1, ictmolecule

				do m = 1, 29

					n = (j - 1) * 29 + m
					rxs(n, ib1000) = rcrdrbufd(ld + 1)
					rys(n, ib1000) = rcrdrbufd(ld + 2)
					rzs(n, ib1000) = rcrdrbufd(ld + 3)
					ld = ld + 3

				enddo

			enddo

		enddo

		do ib1000 = 801, 900

			kd = kd + 1
			icts(ib1000) = ictrbufd(kd)
			ictmolecule = icts(ib1000)

			do j = 1, ictmolecule

				do m = 1, 29

					n = (j - 1) * 29 + m
					rxs(n, ib1000) = rcrdrbufd(ld + 1)
					rys(n, ib1000) = rcrdrbufd(ld + 2)
					rzs(n, ib1000) = rcrdrbufd(ld + 3)
					ld = ld + 3

				enddo

			enddo

		enddo

		do ib512 = 449, 512

			kd = kd + 1
			ictw2(ib512) = ictrbufd(kd)
			ictmolecule = ictw2(ib512)

			do j = 1, ictmolecule

				do m = 1, 3

					n = (j - 1) * 3 + m
					rxw2(n, ib512) = rcrdrbufd(ld + 1)
					ryw2(n, ib512) = rcrdrbufd(ld + 2)
					rzw2(n, ib512) = rcrdrbufd(ld + 3)
					ld = ld + 3

				enddo

			enddo

		enddo

		do ib1000 = 701, 800 

			kd = kd + 1
			icts(ib1000) = ictrbufd(kd)
			ictmolecule = icts(ib1000)

			do j = 1, ictmolecule

				do m = 1, 29

					n = (j - 1) * 29 + m
					rxs(n, ib1000) = rcrdrbufd(ld + 1)
					rys(n, ib1000) = rcrdrbufd(ld + 2)
					rzs(n, ib1000) = rcrdrbufd(ld + 3)
					ld = ld + 3

				enddo

			enddo

		enddo

		do ib512 = 385, 448

			kd = kd + 1
			ictw2(ib512) = ictrbufd(kd)
			ictmolecule = ictw2(ib512)

			do j = 1, ictmolecule

				do m = 1, 3

					n = (j - 1) * 3 + m
					rxw2(n, ib512) = rcrdrbufd(ld + 1)
					ryw2(n, ib512) = rcrdrbufd(ld + 2)
					rzw2(n, ib512) = rcrdrbufd(ld + 3)
					ld = ld + 3

				enddo

			enddo

		enddo

		do ib216 = 181, 216

			kd = kd + 1
			ictw1(ib216) = ictrbufd(kd)
			ictmolecule = ictw1(ib216)

			do j = 1, ictmolecule

				do m = 1, 3

					n = (j - 1) * 3 + m
					rxw1(n, ib216) = rcrdrbufd(ld + 1)
					ryw1(n, ib216) = rcrdrbufd(ld + 2)
					rzw1(n, ib216) = rcrdrbufd(ld + 3)
					ld = ld + 3

				enddo

			enddo

		enddo

		if(ld.ne.ictcrdd) write(*,*) "downSR ld"
		if(kd.ne.464) write(*,*) "downSR kd unload"
		
	endif

	return
	end subroutine
