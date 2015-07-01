	subroutine eval_subtract_fma(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, ib4, ib1000, ib64to1000, ictk, ictn
	integer*4 i, j, k, m, n
	real*8 xm, ym, zm, xmn, ymn, zmn, rmnsq, rmn, fel
	
	do ib4 = 1, 64
	
		ib1000 = ib64to1000(ib4)
		
		ictk = icts(ib1000)
		do i = 1, ictk
		
			do j = 1, 29
			
				m = (i - 1) * 29 + j
				
				xm = rxs(m, ib1000)
				ym = rys(m, ib1000)
				zm = rzs(m, ib1000)
				
				ictn = 1
				k = isub(ictn, j, i, ib4)
				
				do while(k.ne.0)
				
					n = (i - 1) * 29 + k
					
					xmn = xm - rxs(n, ib1000)
					ymn = ym - rys(n, ib1000)
					zmn = zm - rzs(n, ib1000)
					
					rmnsq = xmn**2 + ymn**2 + zmn**2
					rmn = sqrt(rmnsq)
					
					fel = rqs(j) * rqs(k) * xpi / rmn**3
					
					fxs(m, ib1000) = fxs(m, ib1000) - xmn * fel
					fys(m, ib1000) = fys(m, ib1000) - ymn * fel
					fzs(m, ib1000) = fzs(m, ib1000) - zmn * fel
					!write(*,*) "subtracted"
					
					ictn = ictn + 1
					k = isub(ictn, j, i, ib4)
				
				end do
				
			enddo
			
		enddo
	
	enddo
	
	return
	end subroutine
