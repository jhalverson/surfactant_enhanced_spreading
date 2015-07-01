	subroutine nh_water1(raq)
	
	implicit none
	include "mpif.h"
	include "mybox.cmns"
	
	integer*4 ib4, ib216, ib64to216
	integer*4 i, ictk
	
	real*8 raq, sumvelO, sumvelH
	
	sumvelO = 0.0D0
	sumvelH = 0.0D0
	
	do ib4 = 1, 64
	
		ib216 = ib64to216(ib4)
		
		ictk = 3 * ictw1(ib216)
		do i = 1, ictk - 2, 3
		
			sumvelO = sumvelO + vxw1(i, ib4)**2 + vyw1(i, ib4)**2 + vzw1(i, ib4)**2
			
		enddo
		
		do i = 2, ictk - 1, 3
		
			sumvelH = sumvelH + vxw1(i, ib4)**2 + vyw1(i, ib4)**2 + vzw1(i, ib4)**2
			
		enddo
		
		do i = 3, ictk - 0, 3
		
			sumvelH = sumvelH + vxw1(i, ib4)**2 + vyw1(i, ib4)**2 + vzw1(i, ib4)**2
			
		enddo
	
	enddo
	
	raq = emw(1) * sumvelO + emw(2) * sumvelH
	
	return
	end subroutine
