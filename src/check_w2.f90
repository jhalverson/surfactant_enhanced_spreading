	subroutine check_w2(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, ib4, ib512, ib64to512, j, ictk, jstart, ibox2, ibox64
	integer*4 ibO, ibG, ibH, ib64O, ib64G, ib64H
	
	real*8 xO, yO, zO
	real*8 xG, yG, zG
	real*8 xH, yH, zH
	
	do ib4 = 1, 64
	
		ib512 = ib64to512(ib4)
		
		ictk = ictw2(ib512)
		do j = 1, ictk
		
			jstart = (j - 1) * 3
			
			xO = rxw2(jstart + 1, ib512)
			yO = ryw2(jstart + 1, ib512)
			zO = rzw2(jstart + 1, ib512)
			
			xG = rxw2(jstart + 2, ib512)
			yG = ryw2(jstart + 2, ib512)
			zG = rzw2(jstart + 2, ib512)
			
			xH = rxw2(jstart + 3, ib512)
			yH = ryw2(jstart + 3, ib512)
			zH = rzw2(jstart + 3, ib512)
			
			ibO = ibox2(xO, yO, zO)
			ibG = ibox2(xG, yG, zG)
			ibH = ibox2(xH, yH, zH)
			
			if(ibO.eq.ibG.and.ibO.eq.ibH) then
			
				ib64O = ibox64(xO, yO, zO)
				ib64G = ibox64(xG, yG, zG)
				ib64H = ibox64(xH, yH, zH)
				
				if(ib64O.eq.ib64G.and.ib64O.eq.ib64H) then
				
					write(*,*) "not a w2"
				
				endif
				
			endif
			
		enddo
	
	enddo
	
	!if(iaq.eq.0) write(*,*) "check_w2"
	
	return
	end subroutine
