	subroutine init_temp(iaq)
	
	implicit none
	include "mpif.h"
	include "mybox.cmns"
	
	integer*4 iaq, ib4, ib216, ib64to216, ib512, ib64to512, ib1000, ib64to1000
	integer*4 i, ictk, ierr
	integer*4 jctH, jctO, jctSI, jctCH2, jctCH3
	
	real*8 sumvelO, sumvelH, sumvelCH2, sumvelSI, sumvelCH3
	real*8 v2sum, v2total, temp, tscale
	real*8 sumxH, sumyH, sumzH
	real*8 sumxO, sumyO, sumzO
	real*8 sumxSI, sumySI, sumzSI
	real*8 sumxCH2, sumyCH2, sumzCH2
	real*8 sumxCH3, sumyCH3, sumzCH3
	real*8 rnH, rnO, rnSI, rnCH2, rnCH3
	real*8 seed, x, xx, yy, zz, xyz
	
	seed = 0.378389247D0
	x = seed
	
	jctH = 0
	jctO = 0
	jctSI = 0
	jctCH2 = 0
	jctCH3 = 0
	
	sumvelO = 0.0D0
	sumvelH = 0.0D0
	sumvelSI = 0.0D0
	sumvelCH2 = 0.0D0
	sumvelCH3 = 0.0D0
	
	sumxH = 0.0D0
	sumyH = 0.0D0
	sumzH = 0.0D0
	
	sumxO = 0.0D0
	sumyO = 0.0D0
	sumzO = 0.0D0
	
	sumxSI = 0.0D0
	sumySI = 0.0D0
	sumzSI = 0.0D0
	
	sumxCH2 = 0.0D0
	sumyCH2 = 0.0D0
	sumzCH2 = 0.0D0
	
	sumxCH3 = 0.0D0
	sumyCH3 = 0.0D0
	sumzCH3 = 0.0D0
	
	do ib4 = 1, 64
	
		ib216 = ib64to216(ib4)
		ib512 = ib64to512(ib4)
		ib1000 = ib64to1000(ib4)
		
		ictk = 3 * ictw1(ib216)
		do i = 1, ictk
		
			call random(x,xx)
			x = xx
			xx = 2.0D0 * (xx - 0.5D0)
			call random(x, yy)
			x = yy
			yy = 2.0D0 * (yy - 0.5D0)
			call random(x, zz)
			x = zz
			zz = 2.0D0 * (zz - 0.5D0)
			xyz = 1.0D0 / dsqrt(xx * xx + yy * yy + zz * zz)
	   
			vxw1(i, ib4) = xx * xyz
			vyw1(i, ib4) = yy * xyz
			vzw1(i, ib4) = zz * xyz
			
		enddo
		
		ictk = 3 * ictw2(ib512)
		do i = 1, ictk
		
			call random(x,xx)
			x = xx
			xx = 2.0D0 * (xx - 0.5D0)
			call random(x, yy)
			x = yy
			yy = 2.0D0 * (yy - 0.5D0)
			call random(x, zz)
			x = zz
			zz = 2.0D0 * (zz - 0.5D0)
			xyz = 1.0D0 / dsqrt(xx * xx + yy * yy + zz * zz)
	   
			vxw2(i, ib4) = xx * xyz
			vyw2(i, ib4) = yy * xyz
			vzw2(i, ib4) = zz * xyz
			
		enddo
		
		ictk = 29 * icts(ib1000)
		do i = 1, ictk
		
			call random(x,xx)
			x = xx
			xx = 2.0D0 * (xx - 0.5D0)
			call random(x, yy)
			x = yy
			yy = 2.0D0 * (yy - 0.5D0)
			call random(x, zz)
			x = zz
			zz = 2.0D0 * (zz - 0.5D0)
			xyz = 1.0D0 / dsqrt(xx * xx + yy * yy + zz * zz)
			
			vxs(i, ib4) = xx * xyz
			vys(i, ib4) = yy * xyz
			vzs(i, ib4) = zz * xyz
			
		enddo
		
	enddo
	
	do ib4 = 1, 64
	
		ib216 = ib64to216(ib4)
		ib512 = ib64to512(ib4)
		ib1000 = ib64to1000(ib4)
		
		ictk = 3 * ictw1(ib216)
		do i = 1, ictk - 2, 3
		
			sumxO = sumxO + vxw1(i, ib4)
			sumyO = sumyO + vyw1(i, ib4)
			sumzO = sumzO + vzw1(i, ib4)
			jctO = jctO + 1
			
		enddo
		
		do i = 2, ictk - 1, 3
			
			sumxH = sumxH + vxw1(i, ib4)
			sumyH = sumyH + vyw1(i, ib4)
			sumzH = sumzH + vzw1(i, ib4)
			jctH = jctH + 1
			
		enddo
		
		do i = 3, ictk - 0, 3
			
			sumxH = sumxH + vxw1(i, ib4)
			sumyH = sumyH + vyw1(i, ib4)
			sumzH = sumzH + vzw1(i, ib4)
			jctH = jctH + 1
			
		enddo
		
		ictk = 3 * ictw2(ib512)
		do i = 1, ictk - 2, 3
			
			sumxO = sumxO + vxw2(i, ib4)
			sumyO = sumyO + vyw2(i, ib4)
			sumzO = sumzO + vzw2(i, ib4)
			jctO = jctO + 1
			
		enddo
		
		do i = 2, ictk - 1, 3
		
			sumxH = sumxH + vxw2(i, ib4)
			sumyH = sumyH + vyw2(i, ib4)
			sumzH = sumzH + vzw2(i, ib4)
			jctH = jctH + 1
			
		enddo
		
		do i = 3, ictk - 0, 3
		
			sumxH = sumxH + vxw2(i, ib4)
			sumyH = sumyH + vyw2(i, ib4)
			sumzH = sumzH + vzw2(i, ib4)
			jctH = jctH + 1
			
		enddo
		
		ictk = 29 * icts(ib1000)
		do i = 1, ictk - 28, 29
		
			sumxH = sumxH + vxs(i, ib4)
			sumyH = sumyH + vys(i, ib4)
			sumzH = sumzH + vzs(i, ib4)
			jctH = jctH + 1
			
		enddo
		
		do i = 2, ictk - 27, 29
			
			sumxO = sumxO + vxs(i, ib4)
			sumyO = sumyO + vys(i, ib4)
			sumzO = sumzO + vzs(i, ib4)
			jctO = jctO + 1
			
		enddo
		
		do i = 3, ictk - 26, 29
		
			sumxCH2 = sumxCH2 + vxs(i, ib4)
			sumyCH2 = sumyCH2 + vys(i, ib4)
			sumzCH2 = sumzCH2 + vzs(i, ib4)
			jctCH2 = jctCH2 + 1
			
		enddo
		
		do i = 4, ictk - 25, 29
		
			sumxCH2 = sumxCH2 + vxs(i, ib4)
			sumyCH2 = sumyCH2 + vys(i, ib4)
			sumzCH2 = sumzCH2 + vzs(i, ib4)
			jctCH2 = jctCH2 + 1
			
		enddo
		
		do i = 5, ictk - 24, 29
		
			sumxO = sumxO + vxs(i, ib4)
			sumyO = sumyO + vys(i, ib4)
			sumzO = sumzO + vzs(i, ib4)
			jctO = jctO + 1
			
		enddo
		
		do i = 6, ictk - 23, 29
		
			sumxCH2 = sumxCH2 + vxs(i, ib4)
			sumyCH2 = sumyCH2 + vys(i, ib4)
			sumzCH2 = sumzCH2 + vzs(i, ib4)
			jctCH2 = jctCH2 + 1
			
		enddo
		
		do i = 7, ictk - 22, 29
		
			sumxCH2 = sumxCH2 + vxs(i, ib4)
			sumyCH2 = sumyCH2 + vys(i, ib4)
			sumzCH2 = sumzCH2 + vzs(i, ib4)
			jctCH2 = jctCH2 + 1
			
		enddo
		
		do i = 8, ictk - 21, 29
		
			sumxO = sumxO + vxs(i, ib4)
			sumyO = sumyO + vys(i, ib4)
			sumzO = sumzO + vzs(i, ib4)
			jctO = jctO + 1
			
		enddo
		
		do i = 9, ictk - 20, 29
		
			sumxCH2 = sumxCH2 + vxs(i, ib4)
			sumyCH2 = sumyCH2 + vys(i, ib4)
			sumzCH2 = sumzCH2 + vzs(i, ib4)
			jctCH2 = jctCH2 + 1
			
		enddo
		
		do i = 10, ictk - 19, 29
		
			sumxCH2 = sumxCH2 + vxs(i, ib4)
			sumyCH2 = sumyCH2 + vys(i, ib4)
			sumzCH2 = sumzCH2 + vzs(i, ib4)
			jctCH2 = jctCH2 + 1
			
		enddo
		
		do i = 11, ictk - 18, 29
		
			sumxO = sumxO + vxs(i, ib4)
			sumyO = sumyO + vys(i, ib4)
			sumzO = sumzO + vzs(i, ib4)
			jctO = jctO + 1
			
		enddo
		
		do i = 12, ictk - 17, 29
		
			sumxCH2 = sumxCH2 + vxs(i, ib4)
			sumyCH2 = sumyCH2 + vys(i, ib4)
			sumzCH2 = sumzCH2 + vzs(i, ib4)
			jctCH2 = jctCH2 + 1
			
		enddo
		
		do i = 13, ictk - 16, 29
		
			sumxCH2 = sumxCH2 + vxs(i, ib4)
			sumyCH2 = sumyCH2 + vys(i, ib4)
			sumzCH2 = sumzCH2 + vzs(i, ib4)
			jctCH2 = jctCH2 + 1
			
		enddo
		
		do i = 14, ictk - 15, 29
		
			sumxO = sumxO + vxs(i, ib4)
			sumyO = sumyO + vys(i, ib4)
			sumzO = sumzO + vzs(i, ib4)
			jctO = jctO + 1
			
		enddo
		
		do i = 15, ictk - 14, 29
		
			sumxCH2 = sumxCH2 + vxs(i, ib4)
			sumyCH2 = sumyCH2 + vys(i, ib4)
			sumzCH2 = sumzCH2 + vzs(i, ib4)
			jctCH2 = jctCH2 + 1
			
		enddo
		
		do i = 16, ictk - 13, 29
		
			sumxCH2 = sumxCH2 + vxs(i, ib4)
			sumyCH2 = sumyCH2 + vys(i, ib4)
			sumzCH2 = sumzCH2 + vzs(i, ib4)
			jctCH2 = jctCH2 + 1
			
		enddo
		
		do i = 17, ictk - 12, 29
		
			sumxCH2 = sumxCH2 + vxs(i, ib4)
			sumyCH2 = sumyCH2 + vys(i, ib4)
			sumzCH2 = sumzCH2 + vzs(i, ib4)
			jctCH2 = jctCH2 + 1
			
		enddo
		
		do i = 18, ictk - 11, 29
		
			sumxSI = sumxSI + vxs(i, ib4)
			sumySI = sumySI + vys(i, ib4)
			sumzSI = sumzSI + vzs(i, ib4)
			jctSI = jctSI + 1
			
		enddo
		
		do i = 19, ictk - 10, 29
		
			sumxCH3 = sumxCH3 + vxs(i, ib4)
			sumyCH3 = sumyCH3 + vys(i, ib4)
			sumzCH3 = sumzCH3 + vzs(i, ib4)
			jctCH3 = jctCH3 + 1
			
		enddo
		
		do i = 20, ictk - 9, 29
		
			sumxO = sumxO + vxs(i, ib4)
			sumyO = sumyO + vys(i, ib4)
			sumzO = sumzO + vzs(i, ib4)
			jctO = jctO + 1
			
		enddo
		
		do i = 21, ictk - 8, 29
		
			sumxSI = sumxSI + vxs(i, ib4)
			sumySI = sumySI + vys(i, ib4)
			sumzSI = sumzSI + vzs(i, ib4)
			jctSI = jctSI + 1
			
		enddo
		
		do i = 22, ictk - 7, 29
			
			sumxCH3 = sumxCH3 + vxs(i, ib4)
			sumyCH3 = sumyCH3 + vys(i, ib4)
			sumzCH3 = sumzCH3 + vzs(i, ib4)
			jctCH3 = jctCH3 + 1
			
		enddo
		
		do i = 23, ictk - 6, 29
			
			sumxCH3 = sumxCH3 + vxs(i, ib4)
			sumyCH3 = sumyCH3 + vys(i, ib4)
			sumzCH3 = sumzCH3 + vzs(i, ib4)
			jctCH3 = jctCH3 + 1
			
		enddo
		
		do i = 24, ictk - 5, 29
		
			sumxCH3 = sumxCH3 + vxs(i, ib4)
			sumyCH3 = sumyCH3 + vys(i, ib4)
			sumzCH3 = sumzCH3 + vzs(i, ib4)
			jctCH3 = jctCH3 + 1
			
		enddo
		
		do i = 25, ictk - 4, 29
		
			sumxO = sumxO + vxs(i, ib4)
			sumyO = sumyO + vys(i, ib4)
			sumzO = sumzO + vzs(i, ib4)
			jctO = jctO + 1
			
		enddo
		
		do i = 26, ictk - 3, 29
		
			sumxSI = sumxSI + vxs(i, ib4)
			sumySI = sumySI + vys(i, ib4)
			sumzSI = sumzSI + vzs(i, ib4)
			jctSI = jctSI + 1
			
		enddo
		
		do i = 27, ictk - 2, 29
		
			sumxCH3 = sumxCH3 + vxs(i, ib4)
			sumyCH3 = sumyCH3 + vys(i, ib4)
			sumzCH3 = sumzCH3 + vzs(i, ib4)
			jctCH3 = jctCH3 + 1
			
		enddo
		
		do i = 28, ictk - 1, 29
		
			sumxCH3 = sumxCH3 + vxs(i, ib4)
			sumyCH3 = sumyCH3 + vys(i, ib4)
			sumzCH3 = sumzCH3 + vzs(i, ib4)
			jctCH3 = jctCH3 + 1
			
		enddo
		
		do i = 29, ictk - 0, 29
		
			sumxCH3 = sumxCH3 + vxs(i, ib4)
			sumyCH3 = sumyCH3 + vys(i, ib4)
			sumzCH3 = sumzCH3 + vzs(i, ib4)
			jctCH3 = jctCH3 + 1
			
		enddo
	
	enddo
	
	rnH = 1.0D0 / jctH
	rnO = 1.0D0 / jctO
	rnSI = 1.0D0 / jctSI
	rnCH2 = 1.0D0 / jctCH2
	rnCH3 = 1.0D0 / jctCH3
	
	do ib4 = 1, 64
	
		ib216 = ib64to216(ib4)
		ib512 = ib64to512(ib4)
		ib1000 = ib64to1000(ib4)
		
		ictk = 3 * ictw1(ib216)
		do i = 1, ictk - 2, 3
		
			vxw1(i, ib4) = vxw1(i, ib4) - sumxO * rnO
			vyw1(i, ib4) = vyw1(i, ib4) - sumyO * rnO
			vzw1(i, ib4) = vzw1(i, ib4) - sumzO * rnO
			sumvelO = sumvelO + vxw1(i, ib4)**2 + vyw1(i, ib4)**2 + vzw1(i, ib4)**2
			
		enddo
		
		do i = 2, ictk - 1, 3
			
			vxw1(i, ib4) = vxw1(i, ib4) - sumxH * rnH
			vyw1(i, ib4) = vyw1(i, ib4) - sumyH * rnH
			vzw1(i, ib4) = vzw1(i, ib4) - sumzH * rnH
			sumvelH = sumvelH + vxw1(i, ib4)**2 + vyw1(i, ib4)**2 + vzw1(i, ib4)**2
			
		enddo
		
		do i = 3, ictk - 0, 3
			
			vxw1(i, ib4) = vxw1(i, ib4) - sumxH * rnH
			vyw1(i, ib4) = vyw1(i, ib4) - sumyH * rnH
			vzw1(i, ib4) = vzw1(i, ib4) - sumzH * rnH
			sumvelH = sumvelH + vxw1(i, ib4)**2 + vyw1(i, ib4)**2 + vzw1(i, ib4)**2
			
		enddo
		
		ictk = 3 * ictw2(ib512)
		do i = 1, ictk - 2, 3
			
			vxw2(i, ib4) = vxw2(i, ib4) - sumxO * rnO
			vyw2(i, ib4) = vyw2(i, ib4) - sumyO * rnO
			vzw2(i, ib4) = vzw2(i, ib4) - sumzO * rnO
			sumvelO = sumvelO + vxw2(i, ib4)**2 + vyw2(i, ib4)**2 + vzw2(i, ib4)**2
			
		enddo
		
		do i = 2, ictk - 1, 3
		
			vxw2(i, ib4) = vxw2(i, ib4) - sumxH * rnH
			vyw2(i, ib4) = vyw2(i, ib4) - sumyH * rnH
			vzw2(i, ib4) = vzw2(i, ib4) - sumzH * rnH
			sumvelH = sumvelH + vxw2(i, ib4)**2 + vyw2(i, ib4)**2 + vzw2(i, ib4)**2
			
		enddo
		
		do i = 3, ictk - 0, 3
		
			vxw2(i, ib4) = vxw2(i, ib4) - sumxH * rnH
			vyw2(i, ib4) = vyw2(i, ib4) - sumyH * rnH
			vzw2(i, ib4) = vzw2(i, ib4) - sumzH * rnH
			sumvelH = sumvelH + vxw2(i, ib4)**2 + vyw2(i, ib4)**2 + vzw2(i, ib4)**2
			
		enddo
		
		ictk = 29 * icts(ib1000)
		do i = 1, ictk - 28, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxH * rnH
			vys(i, ib4) = vys(i, ib4) - sumyH * rnH
			vzs(i, ib4) = vzs(i, ib4) - sumzH * rnH
			sumvelH = sumvelH + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 2, ictk - 27, 29
			
			vxs(i, ib4) = vxs(i, ib4) - sumxO * rnO
			vys(i, ib4) = vys(i, ib4) - sumyO * rnO
			vzs(i, ib4) = vzs(i, ib4) - sumzO * rnO
			sumvelO = sumvelO + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 3, ictk - 26, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH2 * rnCH2
			vys(i, ib4) = vys(i, ib4) - sumyCH2 * rnCH2
			vzs(i, ib4) = vzs(i, ib4) - sumzCH2 * rnCH2
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 4, ictk - 25, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH2 * rnCH2
			vys(i, ib4) = vys(i, ib4) - sumyCH2 * rnCH2
			vzs(i, ib4) = vzs(i, ib4) - sumzCH2 * rnCH2
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 5, ictk - 24, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxO * rnO
			vys(i, ib4) = vys(i, ib4) - sumyO * rnO
			vzs(i, ib4) = vzs(i, ib4) - sumzO * rnO
			sumvelO = sumvelO + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 6, ictk - 23, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH2 * rnCH2
			vys(i, ib4) = vys(i, ib4) - sumyCH2 * rnCH2
			vzs(i, ib4) = vzs(i, ib4) - sumzCH2 * rnCH2
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 7, ictk - 22, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH2 * rnCH2
			vys(i, ib4) = vys(i, ib4) - sumyCH2 * rnCH2
			vzs(i, ib4) = vzs(i, ib4) - sumzCH2 * rnCH2
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 8, ictk - 21, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxO * rnO
			vys(i, ib4) = vys(i, ib4) - sumyO * rnO
			vzs(i, ib4) = vzs(i, ib4) - sumzO * rnO
			sumvelO = sumvelO + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 9, ictk - 20, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH2 * rnCH2
			vys(i, ib4) = vys(i, ib4) - sumyCH2 * rnCH2
			vzs(i, ib4) = vzs(i, ib4) - sumzCH2 * rnCH2
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 10, ictk - 19, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH2 * rnCH2
			vys(i, ib4) = vys(i, ib4) - sumyCH2 * rnCH2
			vzs(i, ib4) = vzs(i, ib4) - sumzCH2 * rnCH2
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 11, ictk - 18, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxO * rnO
			vys(i, ib4) = vys(i, ib4) - sumyO * rnO
			vzs(i, ib4) = vzs(i, ib4) - sumzO * rnO
			sumvelO = sumvelO + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 12, ictk - 17, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH2 * rnCH2
			vys(i, ib4) = vys(i, ib4) - sumyCH2 * rnCH2
			vzs(i, ib4) = vzs(i, ib4) - sumzCH2 * rnCH2
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 13, ictk - 16, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH2 * rnCH2
			vys(i, ib4) = vys(i, ib4) - sumyCH2 * rnCH2
			vzs(i, ib4) = vzs(i, ib4) - sumzCH2 * rnCH2
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 14, ictk - 15, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxO * rnO
			vys(i, ib4) = vys(i, ib4) - sumyO * rnO
			vzs(i, ib4) = vzs(i, ib4) - sumzO * rnO
			sumvelO = sumvelO + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 15, ictk - 14, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH2 * rnCH2
			vys(i, ib4) = vys(i, ib4) - sumyCH2 * rnCH2
			vzs(i, ib4) = vzs(i, ib4) - sumzCH2 * rnCH2
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 16, ictk - 13, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH2 * rnCH2
			vys(i, ib4) = vys(i, ib4) - sumyCH2 * rnCH2
			vzs(i, ib4) = vzs(i, ib4) - sumzCH2 * rnCH2
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 17, ictk - 12, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH2 * rnCH2
			vys(i, ib4) = vys(i, ib4) - sumyCH2 * rnCH2
			vzs(i, ib4) = vzs(i, ib4) - sumzCH2 * rnCH2
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 18, ictk - 11, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxSI * rnSI
			vys(i, ib4) = vys(i, ib4) - sumySI * rnSI
			vzs(i, ib4) = vzs(i, ib4) - sumzSI * rnSI
			sumvelSI = sumvelSI + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 19, ictk - 10, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH3 * rnCH3
			vys(i, ib4) = vys(i, ib4) - sumyCH3 * rnCH3
			vzs(i, ib4) = vzs(i, ib4) - sumzCH3 * rnCH3
			sumvelCH3 = sumvelCH3 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 20, ictk - 9, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxO * rnO
			vys(i, ib4) = vys(i, ib4) - sumyO * rnO
			vzs(i, ib4) = vzs(i, ib4) - sumzO * rnO
			sumvelO = sumvelO + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 21, ictk - 8, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxSI * rnSI
			vys(i, ib4) = vys(i, ib4) - sumySI * rnSI
			vzs(i, ib4) = vzs(i, ib4) - sumzSI * rnSI
			sumvelSI = sumvelSI + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 22, ictk - 7, 29
			
			vxs(i, ib4) = vxs(i, ib4) - sumxCH3 * rnCH3
			vys(i, ib4) = vys(i, ib4) - sumyCH3 * rnCH3
			vzs(i, ib4) = vzs(i, ib4) - sumzCH3 * rnCH3
			sumvelCH3 = sumvelCH3 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 23, ictk - 6, 29
			
			vxs(i, ib4) = vxs(i, ib4) - sumxCH3 * rnCH3
			vys(i, ib4) = vys(i, ib4) - sumyCH3 * rnCH3
			vzs(i, ib4) = vzs(i, ib4) - sumzCH3 * rnCH3
			sumvelCH3 = sumvelCH3 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 24, ictk - 5, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH3 * rnCH3
			vys(i, ib4) = vys(i, ib4) - sumyCH3 * rnCH3
			vzs(i, ib4) = vzs(i, ib4) - sumzCH3 * rnCH3
			sumvelCH3 = sumvelCH3 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 25, ictk - 4, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxO * rnO
			vys(i, ib4) = vys(i, ib4) - sumyO * rnO
			vzs(i, ib4) = vzs(i, ib4) - sumzO * rnO
			sumvelO = sumvelO + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 26, ictk - 3, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxSI * rnSI
			vys(i, ib4) = vys(i, ib4) - sumySI * rnSI
			vzs(i, ib4) = vzs(i, ib4) - sumzSI * rnSI
			sumvelSI = sumvelSI + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 27, ictk - 2, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH3 * rnCH3
			vys(i, ib4) = vys(i, ib4) - sumyCH3 * rnCH3
			vzs(i, ib4) = vzs(i, ib4) - sumzCH3 * rnCH3
			sumvelCH3 = sumvelCH3 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 28, ictk - 1, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH3 * rnCH3
			vys(i, ib4) = vys(i, ib4) - sumyCH3 * rnCH3
			vzs(i, ib4) = vzs(i, ib4) - sumzCH3 * rnCH3
			sumvelCH3 = sumvelCH3 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 29, ictk - 0, 29
		
			vxs(i, ib4) = vxs(i, ib4) - sumxCH3 * rnCH3
			vys(i, ib4) = vys(i, ib4) - sumyCH3 * rnCH3
			vzs(i, ib4) = vzs(i, ib4) - sumzCH3 * rnCH3
			sumvelCH3 = sumvelCH3 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
		
		enddo
	
	enddo
	
	v2sum = emw(1) * sumvelO + emw(2) * sumvelH + ems(3) * sumvelCH2 + ems(18) * sumvelSI + ems(19) * sumvelCH3
	call mpi_allreduce(v2sum, v2total, 1, mpi_real8, mpi_sum, mpi_comm_world, ierr)
	temp = v2total / dtsq / (3 * (3 * iwater + 29 * isurfactant) - (3 * iwater + 28 * isurfactant))
	tscale = sqrt(tr / temp)
	
	do ib4 = 1, 64
	
		ib216 = ib64to216(ib4)
		ib512 = ib64to512(ib4)
		ib1000 = ib64to1000(ib4)
		
		ictk = 3 * ictw1(ib216)
		do i = 1, ictk
		
			vxw1(i, ib4) = tscale * vxw1(i, ib4)
			vyw1(i, ib4) = tscale * vyw1(i, ib4)
			vzw1(i, ib4) = tscale * vzw1(i, ib4) - 0.0005D0*0
		
		enddo
		
		ictk = 3 * ictw2(ib512)
		do i = 1, ictk
		
			vxw2(i, ib4) = tscale * vxw2(i, ib4)
			vyw2(i, ib4) = tscale * vyw2(i, ib4)
			vzw2(i, ib4) = tscale * vzw2(i, ib4) - 0.0005D0*0
		
		enddo
		
		ictk = 29 * icts(ib1000)
		do i = 1, ictk
		
			vxs(i, ib4) = tscale * vxs(i, ib4)
			vys(i, ib4) = tscale * vys(i, ib4)
			vzs(i, ib4) = tscale * vzs(i, ib4) - 0.0005D0*0
		
		enddo
		
	enddo
	
	if(iaq.eq.0) write(*,*) "Velocities follow a Maxell-Boltzmann distribution."
	
	return
	end subroutine
