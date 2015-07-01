	subroutine computeOnm4(iaq)

	implicit none
	include "mpif.h"
	include "mybox.cmns"

	integer*4 iaq, ierr, ib512, ib64to512
	integer*4 ictOnm4, ib1000, ib64to1000, ictk
	integer*4 ib216, ib64to216
	integer*4 i, j, k, l, m, n, jmolecule, jtype
	integer*4 igetIndex, icoeff, icoeff2
	real*8 rrho(500), ralpha(500), rbeta(500), rqq(500)
	real*8 rctr4x, rctr4y, rctr4z
	complex*16 Ylm
	real*8 qcharge, qchargetotal
	integer*4 iogn, ihgn, iogntotal, ihgntotal

	qcharge = 0.0D0
	iogn = 0
	ihgn = 0

	cO4 = cmplx(0.0D0, 0.0D0)

	do i = 1, 64

		ib216 = ib64to216(i)
		ib512 = ib64to512(i)
		ib1000 = ib64to1000(i)

		ictOnm4 = 0

		rctr4x = rcen4x(i)
		rctr4y = rcen4y(i)
		rctr4z = rcen4z(i)

		! W1
		ictk = 3 * ictw1(ib216)
		do j = 1, ictk

			ictOnm4 = ictOnm4 + 1
			rrho(ictOnm4) = sqrt((rxw1(j, ib216) - rctr4x)**2 + (ryw1(j, ib216) - rctr4y)**2 + (rzw1(j, ib216) - rctr4z)**2)
			ralpha(ictOnm4) = acos((rzw1(j, ib216) - rctr4z) / rrho(ictOnm4))
			rbeta(ictOnm4) = atan2(ryw1(j, ib216) - rctr4y, rxw1(j, ib216) - rctr4x)
			jmolecule = (j - 1) / 3
			jtype = j - 3 * jmolecule
			rqq(ictOnm4) = rqw(jtype)
			qcharge = qcharge + rqq(ictOnm4)
			if(jtype.eq.1) then
				iogn = iogn + 1
			else
				ihgn = ihgn + 1
			endif
		enddo

		! W2
		ictk = ictW2Central(i)
		do j = 1, ictk

			ictOnm4 = ictOnm4 + 1
			m = W2Central(j, 3, i)
			rrho(ictOnm4) = sqrt((rxw2(m, ib512) - rctr4x)**2 + (ryw2(m, ib512) - rctr4y)**2 + (rzw2(m, ib512) - rctr4z)**2)
			ralpha(ictOnm4) = acos((rzw2(m, ib512) - rctr4z) / rrho(ictOnm4))
			rbeta(ictOnm4) = atan2(ryw2(m, ib512) - rctr4y, rxw2(m, ib512) - rctr4x)
			jtype = W2Central(j, 4, i)
			rqq(ictOnm4) = rqw(jtype)
			qcharge = qcharge + rqq(ictOnm4)
			if(jtype.eq.1) then
				iogn = iogn + 1
			else
				ihgn = ihgn + 1
			endif

		enddo

		! S
		ictk = ictSCentral(i)
		do j = 1, ictk

			ictOnm4 = ictOnm4 + 1
			m = SCentral(j, 3, i)
			rrho(ictOnm4) = sqrt((rxs(m, ib1000) - rctr4x)**2 + (rys(m, ib1000) - rctr4y)**2 + (rzs(m, ib1000) - rctr4z)**2)
			ralpha(ictOnm4) = acos((rzs(m, ib1000) - rctr4z) / rrho(ictOnm4))
			rbeta(ictOnm4) = atan2(rys(m, ib1000) - rctr4y, rxs(m, ib1000) - rctr4x)
			jtype = SCentral(j, 4, i)
			rqq(ictOnm4) = rqs(jtype)
			qcharge = qcharge + rqq(ictOnm4)

		enddo

		! W2n
		ictk = ictW2nCentral(i)
		do j = 1, ictk

			ictOnm4 = ictOnm4 + 1
			l = W2nCentral(j, 1, i)
			m = W2nCentral(j, 3, i)
			rrho(ictOnm4) = sqrt((rxw2(m, l) - rctr4x)**2 + (ryw2(m, l) - rctr4y)**2 + (rzw2(m, l) - rctr4z)**2)
			ralpha(ictOnm4) = acos((rzw2(m, l) - rctr4z) / rrho(ictOnm4))
			rbeta(ictOnm4) = atan2(ryw2(m, l) - rctr4y, rxw2(m, l) - rctr4x)
			jtype = W2nCentral(j, 4, i)
			if(jtype.eq.1) write(*,*) "computeOnm4 jtype"
			rqq(ictOnm4) = rqw(jtype)
			qcharge = qcharge + rqq(ictOnm4)
			ihgn = ihgn + 1

		enddo

		! Sn
		ictk = ictSnCentral(i)
		do j = 1, ictk

			ictOnm4 = ictOnm4 + 1
			l = SnCentral(j, 1, i)
			m = SnCentral(j, 3, i)
			rrho(ictOnm4) = sqrt((rxs(m, l) - rctr4x)**2 + (rys(m, l) - rctr4y)**2 + (rzs(m, l) - rctr4z)**2)
			ralpha(ictOnm4) = acos((rzs(m, l) - rctr4z) / rrho(ictOnm4))
			rbeta(ictOnm4) = atan2(rys(m, l) - rctr4y, rxs(m, l) - rctr4x)
			jtype = SnCentral(j, 4, i)
			rqq(ictOnm4) = rqs(jtype)
			qcharge = qcharge + rqq(ictOnm4)

		enddo

		! Snn
		ictk = ictSnnCentral(i)
		do j = 1, ictk

			ictOnm4 = ictOnm4 + 1
			l = SnnCentral(j, 1, i)
			m = SnnCentral(j, 3, i)
			rrho(ictOnm4) = sqrt((rxs(m, l) - rctr4x)**2 + (rys(m, l) - rctr4y)**2 + (rzs(m, l) - rctr4z)**2)
			ralpha(ictOnm4) = acos((rzs(m, l) - rctr4z) / rrho(ictOnm4))
			rbeta(ictOnm4) = atan2(rys(m, l) - rctr4y, rxs(m, l) - rctr4x)
			jtype = SnnCentral(j, 4, i)
			rqq(ictOnm4) = rqs(jtype)
			qcharge = qcharge + rqq(ictOnm4)

		enddo

		!if(iaq.eq.13) write(*,*) "charge = ", qcharge, "ib4 = ", i

		do n = 0, p

			do m = 0, n

				icoeff = igetIndex(n, m)

				do j = 1, ictOnm4

					cO4(icoeff, ib512) = cO4(icoeff, ib512) + rqq(j) * rrho(j)**n * Ylm(n, -m, ralpha(j), rbeta(j))

				enddo

				if(m.gt.0) then

					icoeff2 = igetIndex(n, -m)
					cO4(icoeff2, ib512) = conjg(cO4(icoeff, ib512))
				endif

			enddo

		enddo
		
		if(.false.) then
		
		if(iaq.eq.13.and.i.le.4) then
		write(*,*) "******", i
		write(*,'(f8.3,f9.3,f9.3)') rctr4x, rctr4y, rctr4z
		do m = 1, lp
		write(*,'(f15.3,f15.3)') cO4(m, ib512)
		enddo
		endif

		if(iaq.eq.13.and.i.ge.5.and.i.le.8) then
                write(*,*) "******", i
                write(*,'(f8.3,f9.3,f9.3)') rctr4x, rctr4y, rctr4z
                do m = 1, lp
                write(*,'(f15.3,f15.3)') cO4(m, ib512)
                enddo
                endif

		if(iaq.eq.13.and.i.ge.17.and.i.le.20) then
                write(*,*) "******", i
                write(*,'(f8.3,f9.3,f9.3)') rctr4x, rctr4y, rctr4z
                do m = 1, lp
                write(*,'(f15.3,f15.3)') cO4(m, ib512)
                enddo
                endif

		if(iaq.eq.13.and.i.ge.21.and.i.le.24) then
                write(*,*) "******", i
                write(*,'(f8.3,f9.3,f9.3)') rctr4x, rctr4y, rctr4z
                do m = 1, lp
                write(*,'(f15.3,f15.3)') cO4(m, ib512)
                enddo
                endif

		endif

	enddo

	!write(*,*) "mype=",iaq,"q=",qcharge
	call mpi_reduce(qcharge, qchargetotal, 1, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
	!if(iaq.eq.0) write(*,*) "total charge=", qchargetotal
	call mpi_reduce(iogn, iogntotal, 1, mpi_integer4, mpi_sum, 0, mpi_comm_world, ierr)
	!if(iaq.eq.0) write(*,*) "total oxygen=", iogntotal
	call mpi_reduce(ihgn, ihgntotal, 1, mpi_integer4, mpi_sum, 0, mpi_comm_world, ierr)
	!if(iaq.eq.0) write(*,*) "total hydrogen=", ihgntotal

	if(iaq.eq.0.and.qchargetotal**2.gt.0.001D0) write(*,*) "computeOnm4 q-squared", qchargetotal
	if(iaq.eq.0.and.ihgntotal / iogntotal.ne.2) write(*,*) "computeOnm4 O:H ratio", iogntotal, ihgntotal
	if(iaq.eq.0.and.iogntotal.ne.iwater) write(*,*) "computeOnm4 missing_oxygen", iogntotal, iwater

	return
	end subroutine
