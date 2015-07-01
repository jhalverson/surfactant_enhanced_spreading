	subroutine evalLjk4(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 comm3Dq, iaq
	integer*4 ib4, ib1000, ib64to1000, ictk
	integer*4 ib512, ib64to512, ib216, ib64to216
	integer*4 i, j, k, m, n, mn, icoeff, igetIndex
	real*8 rctr4x, rctr4y, rctr4z
	real*8 rxpb, rypb, rzpb, rr, rtheta, rphi
	real*8 rp1, rp3, rq1, rq2, rq3, rq5
	real*8 rfctr, rLegendreP, rqq
	complex*16 Ylm, zn, rp2, rq0, rq4, rq6

	zn = cmplx(0.0D0, 1.0D0)

	!if(iaq.eq.0) write(*,*) "return from Ljk4"
	!return
	
	do ib4 = 1, 64

		ib216 = ib64to216(ib4)
		ib512 = ib64to512(ib4)
		ib1000 = ib64to1000(ib4)

		rctr4x = rcen4x(ib4)
		rctr4y = rcen4y(ib4)
		rctr4z = rcen4z(ib4)

		! W1
		ictk = 3 * ictw1(ib216)
		do i = 1, ictk

			rqq = rqw(2) * xpi
			if(mod(i + 2, 3).eq.0) rqq = rqw(1) * xpi

			rxpb = rxw1(i, ib216) - rctr4x
			rypb = ryw1(i, ib216) - rctr4y
			rzpb = rzw1(i, ib216) - rctr4z

			rr = sqrt(rxpb**2 + rypb**2 + rzpb**2)
			rtheta = acos(rzpb / rr)
			rphi = atan2(rypb, rxpb)
			rp1 = rzpb / rr
			rp2 = (rxpb + rypb * zn) / sqrt(rxpb**2 + rypb**2)
			rp3 = 1.0D0 / (rxpb**2 + rypb**2)

			do j = 0, p

				do k = -j, j

					icoeff = igetIndex(j, k)
					rq0 = cL4(icoeff, ib4) * rfctr(j, k) * rp2**k * rr**(j - 1)
					rq1 = (j + abs(k)) * rLegendreP(j - 1, abs(k), rp1)
					rq2 = rLegendreP(j, abs(k), rp1)
					rq3 = rxpb * rzpb * rq1
					rq4 = (zn * k * rypb - j * rxpb) * rr * rq2
					rq5 = rypb * rzpb * rq1
					rq6 = -(zn * k * rxpb + j * rypb) * rr * rq2

					fxw1(i, ib216) = fxw1(i, ib216) + rqq * rp3 * rq0 * (rq3 + rq4)
					fyw1(i, ib216) = fyw1(i, ib216) + rqq * rp3 * rq0 * (rq5 + rq6)
					fzw1(i, ib216) = fzw1(i, ib216) - rqq * rq0 * rq1

				enddo

			enddo

		enddo

		! w2
		ictk = ictW2Central(ib4)
		do i = 1, ictk

			m = W2Central(i, 3, ib4)
			n = W2Central(i, 4, ib4)
			rqq = rqw(n) * xpi

			rxpb = rxw2(m, ib512) - rctr4x
			rypb = ryw2(m, ib512) - rctr4y
			rzpb = rzw2(m, ib512) - rctr4z

			rr = sqrt(rxpb**2 + rypb**2 + rzpb**2)
			rtheta = acos(rzpb / rr)
			rphi = atan2(rypb, rxpb)
			rp1 = rzpb / rr
			rp2 = (rxpb + rypb * zn) / sqrt(rxpb**2 + rypb**2)
			rp3 = 1.0D0 / (rxpb**2 + rypb**2)

			do j = 0, p

				do k = -j, j

					icoeff = igetIndex(j, k)
					rq0 = cL4(icoeff, ib4) * rfctr(j, k) * rp2**k * rr**(j - 1)
					rq1 = (j + abs(k)) * rLegendreP(j - 1, abs(k), rp1)
					rq2 = rLegendreP(j, abs(k), rp1)
					rq3 = rxpb * rzpb * rq1
					rq4 = (zn * k * rypb - j * rxpb) * rr * rq2
					rq5 = rypb * rzpb * rq1
					rq6 = -(zn * k * rxpb + j * rypb) * rr * rq2

					fxw2(m, ib512) = fxw2(m, ib512) + rqq * rp3 * rq0 * (rq3 + rq4)
					fyw2(m, ib512) = fyw2(m, ib512) + rqq * rp3 * rq0 * (rq5 + rq6)
					fzw2(m, ib512) = fzw2(m, ib512) - rqq * rq0 * rq1

				enddo

			enddo

		enddo
		
		! W2n
		ictk = ictW2nCentral(ib4)
		do i = 1, ictk

			mn = W2nCentral(i, 1, ib4)
			m = W2nCentral(i, 3, ib4)
			rqq = rqw(2) * xpi

			rxpb = rxw2(m, mn) - rctr4x
			rypb = ryw2(m, mn) - rctr4y
			rzpb = rzw2(m, mn) - rctr4z

			rr = sqrt(rxpb**2 + rypb**2 + rzpb**2)
			rtheta = acos(rzpb / rr)
			rphi = atan2(rypb, rxpb)
			rp1 = rzpb / rr
			rp2 = (rxpb + rypb * zn) / sqrt(rxpb**2 + rypb**2)
			rp3 = 1.0D0 / (rxpb**2 + rypb**2)

			do j = 0, p

				do k = -j, j

					icoeff = igetIndex(j, k)
					rq0 = cL4(icoeff, ib4) * rfctr(j, k) * rp2**k * rr**(j - 1)
					rq1 = (j + abs(k)) * rLegendreP(j - 1, abs(k), rp1)
					rq2 = rLegendreP(j, abs(k), rp1)
					rq3 = rxpb * rzpb * rq1
					rq4 = (zn * k * rypb - j * rxpb) * rr * rq2
					rq5 = rypb * rzpb * rq1
					rq6 = -(zn * k * rxpb + j * rypb) * rr * rq2

					fxw2(m, mn) = fxw2(m, mn) + rqq * rp3 * rq0 * (rq3 + rq4)
					fyw2(m, mn) = fyw2(m, mn) + rqq * rp3 * rq0 * (rq5 + rq6)
					fzw2(m, mn) = fzw2(m, mn) - rqq * rq0 * rq1

				enddo

			enddo

		enddo

		! S
		ictk = ictSCentral(ib4)
		do i = 1, ictk

			m = SCentral(i, 3, ib4)
			n = SCentral(i, 4, ib4)
			rqq = rqs(n) * xpi

			rxpb = rxs(m, ib1000) - rctr4x
			rypb = rys(m, ib1000) - rctr4y
			rzpb = rzs(m, ib1000) - rctr4z

			rr = sqrt(rxpb**2 + rypb**2 + rzpb**2)
			rtheta = acos(rzpb / rr)
			rphi = atan2(rypb, rxpb)
			rp1 = rzpb / rr
			rp2 = (rxpb + rypb * zn) / sqrt(rxpb**2 + rypb**2)
			rp3 = 1.0D0 / (rxpb**2 + rypb**2)

			do j = 0, p

				do k = -j, j

					icoeff = igetIndex(j, k)
					rq0 = cL4(icoeff, ib4) * rfctr(j, k) * rp2**k * rr**(j - 1)
					rq1 = (j + abs(k)) * rLegendreP(j - 1, abs(k), rp1)
					rq2 = rLegendreP(j, abs(k), rp1)
					rq3 = rxpb * rzpb * rq1
					rq4 = (zn * k * rypb - j * rxpb) * rr * rq2
					rq5 = rypb * rzpb * rq1
					rq6 = -(zn * k * rxpb + j * rypb) * rr * rq2

					fxs(m, ib1000) = fxs(m, ib1000) + rqq * rp3 * rq0 * (rq3 + rq4)
					fys(m, ib1000) = fys(m, ib1000) + rqq * rp3 * rq0 * (rq5 + rq6)
					fzs(m, ib1000) = fzs(m, ib1000) - rqq * rq0 * rq1

				enddo

			enddo

		enddo
		
		! Sn
		ictk = ictSnCentral(ib4)
		do i = 1, ictk
		
			mn = SnCentral(i, 1, ib4)
			m  = SnCentral(i, 3, ib4)
			n  = SnCentral(i, 4, ib4)
			rqq = rqs(n) * xpi
			
			rxpb = rxs(m, mn) - rctr4x
			rypb = rys(m, mn) - rctr4y
			rzpb = rzs(m, mn) - rctr4z

			rr = sqrt(rxpb**2 + rypb**2 + rzpb**2)
			rtheta = acos(rzpb / rr)
			rphi = atan2(rypb, rxpb)
			rp1 = rzpb / rr
			rp2 = (rxpb + rypb * zn) / sqrt(rxpb**2 + rypb**2)
			rp3 = 1.0D0 / (rxpb**2 + rypb**2)

			do j = 0, p

				do k = -j, j

					icoeff = igetIndex(j, k)
					rq0 = cL4(icoeff, ib4) * rfctr(j, k) * rp2**k * rr**(j - 1)
					rq1 = (j + abs(k)) * rLegendreP(j - 1, abs(k), rp1)
					rq2 = rLegendreP(j, abs(k), rp1)
					rq3 = rxpb * rzpb * rq1
					rq4 = (zn * k * rypb - j * rxpb) * rr * rq2
					rq5 = rypb * rzpb * rq1
					rq6 = -(zn * k * rxpb + j * rypb) * rr * rq2

					fxs(m, mn) = fxs(m, mn) + rqq * rp3 * rq0 * (rq3 + rq4)
					fys(m, mn) = fys(m, mn) + rqq * rp3 * rq0 * (rq5 + rq6)
					fzs(m, mn) = fzs(m, mn) - rqq * rq0 * rq1

				enddo

			enddo

		enddo

		! Snn
		ictk = ictSnnCentral(ib4)
		do i = 1, ictk

			mn = SnnCentral(i, 1, ib4)
			m  = SnnCentral(i, 3, ib4)
			n  = SnnCentral(i, 4, ib4)
			rqq = rqs(n) * xpi

			rxpb = rxs(m, mn) - rctr4x
			rypb = rys(m, mn) - rctr4y
			rzpb = rzs(m, mn) - rctr4z

			rr = sqrt(rxpb**2 + rypb**2 + rzpb**2)
			rtheta = acos(rzpb / rr)
			rphi = atan2(rypb, rxpb)
			rp1 = rzpb / rr
			rp2 = (rxpb + rypb * zn) / sqrt(rxpb**2 + rypb**2)
			rp3 = 1.0D0 / (rxpb**2 + rypb**2)

			do j = 0, p

				do k = -j, j

					icoeff = igetIndex(j, k)
					rq0 = cL4(icoeff, ib4) * rfctr(j, k) * rp2**k * rr**(j - 1)
					rq1 = (j + abs(k)) * rLegendreP(j - 1, abs(k), rp1)
					rq2 = rLegendreP(j, abs(k), rp1)
					rq3 = rxpb * rzpb * rq1
					rq4 = (zn * k * rypb - j * rxpb) * rr * rq2
					rq5 = rypb * rzpb * rq1
					rq6 = -(zn * k * rxpb + j * rypb) * rr * rq2

					fxs(m, mn) = fxs(m, mn) + rqq * rp3 * rq0 * (rq3 + rq4)
					fys(m, mn) = fys(m, mn) + rqq * rp3 * rq0 * (rq5 + rq6)
					fzs(m, mn) = fzs(m, mn) - rqq * rq0 * rq1

				enddo

			enddo

		enddo

	enddo

	return
	end subroutine
