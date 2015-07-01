	subroutine nh_water2(raq)
	
	implicit none
	include "mpif.h"
	include "mybox.cmns"
	
	integer*4 ib4, ib512, ib64to512
	integer*4 i, ictk
	
	real*8 raq, sumvelO, sumvelH
	
	sumvelO = 0.0D0
	sumvelH = 0.0D0
	
	do ib4 = 1, 64
	
		ib512 = ib64to512(ib4)
		
		ictk = 3 * ictw2(ib512)
		do i = 1, ictk - 2, 3
		
			sumvelO = sumvelO + vxw2(i, ib4)**2 + vyw2(i, ib4)**2 + vzw2(i, ib4)**2
			
		enddo
		
		do i = 2, ictk - 1, 3
		
			sumvelH = sumvelH + vxw2(i, ib4)**2 + vyw2(i, ib4)**2 + vzw2(i, ib4)**2
			
		enddo
		
		do i = 3, ictk - 0, 3
		
			sumvelH = sumvelH + vxw2(i, ib4)**2 + vyw2(i, ib4)**2 + vzw2(i, ib4)**2
			
		enddo
	
	enddo
	
	raq = emw(1) * sumvelO + emw(2) * sumvelH
	
	return
	end subroutine
