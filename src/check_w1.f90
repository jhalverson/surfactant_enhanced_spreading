	subroutine check_w1(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, ib4, ib216, ib64to216, j, ictk, jstart, ibox2, ibox64
	integer*4 ibO, ibG, ibH, ib64O, ib64G, ib64H
	
	real*8 xO, yO, zO
	real*8 xG, yG, zG
	real*8 xH, yH, zH
	
	do ib4 = 1, 64
	
		ib216 = ib64to216(ib4)
		
		ictk = ictw1(ib216)
		do j = 1, ictk
		
			jstart = (j - 1) * 3
			
			xO = rxw1(jstart + 1, ib216)
			yO = ryw1(jstart + 1, ib216)
			zO = rzw1(jstart + 1, ib216)
			
			xG = rxw1(jstart + 2, ib216)
			yG = ryw1(jstart + 2, ib216)
			zG = rzw1(jstart + 2, ib216)
			
			xH = rxw1(jstart + 3, ib216)
			yH = ryw1(jstart + 3, ib216)
			zH = rzw1(jstart + 3, ib216)
			
			ibO = ibox2(xO, yO, zO)
			ibG = ibox2(xG, yG, zG)
			ibH = ibox2(xH, yH, zH)
			
			if(ibO.eq.iaq + 1) then
			
				if(ibO.eq.ibG.and.ibO.eq.ibH) then
				
					ib64O = ibox64(xO, yO, zO)
					ib64G = ibox64(xG, yG, zG)
					ib64H = ibox64(xH, yH, zH)
					
					if(ib64O.eq.ib4) then
					
						if(ib64O.eq.ib64G.and.ib64O.eq.ib64H) then
						
						else
						
							write(*,'(i5,i5,i5,i5,i5,i5,i5,i5)') ibO, ibG, ibH, ib64O, ib64G, ib64H, ib4, j
							write(*,*) "not a w1", ictk
						
						endif
						
					else
					
						write(*,'(i5,i5,i5,i5,i5,i5)') ibO, ibG, ibH, ib64O, ib4, j
						write(*,*) "not a w1"
						
					endif
					
				else
				
					write(*,'(i5,i5,i5,i5,i5)') ibO, ibG, ibH, ib4, j
					write(*,*) "not a w1"
					
				endif
				
			else
			
				write(*,'(i5,i5)') ibO, iaq
				write(*,*) "not a w1"
			
			endif
			
		enddo
	
	enddo
	
	!if(iaq.eq.0) write(*,*) "check_w1"
	
	return
	end subroutine
