	subroutine buildinternal(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, ib4, i, j, k, m, n
	integer*4 ib1000, ictk, ibgrand, ictn, ibox1000, ib64to1000
	real*8 xm, ym, zm, rctr1000x, rctr1000y, rctr1000z, xn, yn, zn
	logical bool26, isNearest

	isub = 0
	
	do ib4 = 1, 64
	
		ib1000 = ib64to1000(ib4)
		
		ictk = icts(ib1000)
		do i = 1, ictk
		
			do j = 1, 29
			
				m = (i - 1) * 29 + j
				
				xm = rxs(m, ib1000)
				ym = rys(m, ib1000)
				zm = rzs(m, ib1000)
				
				ibgrand = ibox1000(xm, ym, zm)
				
				rctr1000x = rcen1000x(ibgrand)
				rctr1000y = rcen1000y(ibgrand)
				rctr1000z = rcen1000z(ibgrand)
				
				ictn = 0
				
				do k = 1, 29
				
					n = (i - 1) * 29 + k
					
					if(m.ne.n) then
					
						xn = rxs(n, ib1000)
						yn = rys(n, ib1000)
						zn = rzs(n, ib1000)
						
						bool26 = isNearest(xn, yn, zn, rctr1000x, rctr1000y, rctr1000z)
						
						if(.not.bool26) then
						
							ictn = ictn + 1
							isub(ictn, j, i, ib4) = k
						
						endif
						
					endif
					
				enddo
				
			enddo
			
		enddo
	
	enddo
	
	if(iaq.eq.0) then
	do ib4 = 33, 64
	do i = 1, 1
	do j = 1, 29
	do k = 1, 10
	!write(*,'(i5,i5,i5,i5,i5)') isub(k, j, i, ib4), k, j, i, ib4
	enddo
	enddo
	enddo
	enddo
	endif
	
	return
	end subroutine
