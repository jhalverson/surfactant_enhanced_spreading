	subroutine delete_s(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq
	integer*4 i, j, k, l, g, h
	integer*4 ib4, ib1000, ib64to1000, ibox2, ibox64, ibox64mod
	integer*4 ictk, istart, jstart
	integer*4 ib, ibm, ib64, ib64m, ib1000m
	
	real*8 xOkey, yOkey, zOkey
	real*8 rxL2, ryL2, rzL2
	real*8 rx1000m, ry1000m, rz1000m
	
	real*8 rsurfx(29 * kdeletes), rsurfy(29 * kdeletes), rsurfz(29 * kdeletes)
	real*8 vsurfx(29 * kdeletes), vsurfy(29 * kdeletes), vsurfz(29 * kdeletes)

	jdeleteds = 0
	g = 0
	h = 0
	
	do ib4 = 1, 64
	
		k = 0
		ib1000 = ib64to1000(ib4)
		
		ictk = icts(ib1000)
		do j = 1, ictk
		
			jstart = (j - 1) * 29
			
			xOkey = rxs(jstart + 14, ib1000)
			yOkey = rys(jstart + 14, ib1000)
			zOkey = rzs(jstart + 14, ib1000)
			
			ib = ibox2(xOkey, yOkey, zOkey)
			
			if(ib.eq.iaq + 1) then
			
				ib64 = ibox64(xOkey, yOkey, zOkey)
				
				if(ib64.eq.ib4) then
				
					do l = 1, 29
				
						k = k + 1
						rsurfx(k) = rxs(jstart + l, ib1000)
						rsurfy(k) = rys(jstart + l, ib1000)
						rsurfz(k) = rzs(jstart + l, ib1000)
						
						vsurfx(k) = vxs(jstart + l, ib4)
						vsurfy(k) = vys(jstart + l, ib4)
						vsurfz(k) = vzs(jstart + l, ib4)
						
					enddo
					
				! add to proper internal list
				else
				
					ib1000m = ib64to1000(ib64)
					istart = icts(ib1000m) * 29
					
					do l = 1, 29
					
						rxs(istart + l, ib1000m) = rxs(jstart + l, ib1000)
						rys(istart + l, ib1000m) = rys(jstart + l, ib1000)
						rzs(istart + l, ib1000m) = rzs(jstart + l, ib1000)
						
						vxs(istart + l, ib64) = vxs(jstart + l, ib4)
						vys(istart + l, ib64) = vys(jstart + l, ib4)
						vzs(istart + l, ib64) = vzs(jstart + l, ib4)
						
					enddo
					
					icts(ib1000m) = icts(ib1000m) + 1
					
				endif
			
			! molecule no longer in L2 box, add to ewnsud box
			else
			
				jdeleteds = jdeleteds + 1
			
				rxL2 = rcen2x(ib)
				ryL2 = rcen2y(ib)
				rzL2 = rcen2z(ib)
				
				ib64m = ibox64mod(xOkey, yOkey, zOkey, rxL2, ryL2, rzL2)
				
				isewnsud(g + 1) = ib - 1
				isewnsud(g + 2) = ib64m
				g = g + 2
				
				do l = 1, 29
				
					rsewnsud(h + 1) = rxs(jstart + l, ib1000)
					rsewnsud(h + 2) = rys(jstart + l, ib1000)
					rsewnsud(h + 3) = rzs(jstart + l, ib1000)
					
					rsewnsud(h + 4) = vxs(jstart + l, ib4)
					rsewnsud(h + 5) = vys(jstart + l, ib4)
					rsewnsud(h + 6) = vzs(jstart + l, ib4)
					
					h = h + 6
					
				enddo
			
			endif
			
		enddo
		
		! assign new list
		icts(ib1000) = k / 29
		do i = 1, k
		
			rxs(i, ib1000) = rsurfx(i)
			rys(i, ib1000) = rsurfy(i)
			rzs(i, ib1000) = rsurfz(i)
			
			vxs(i, ib4) = vsurfx(i)
			vys(i, ib4) = vsurfy(i)
			vzs(i, ib4) = vsurfz(i)
			
		enddo

		!if(ictk.ne.k/29) write(*,'(i5,i5,i5,i5,i5)') iaq, ib4, ictk, k / 29, 3

		if(k.gt.29 * kdeletes) write(*,*) "delete_s k"
		if(icts(ib1000).gt.ns) write(*,*) "delete_s icts"
	
	enddo

	!if(jdeleteds.ne.0) write(*,*) "s", iaq, jdeleteds
	if(h.gt.kdeletes296) write(*,*) "delete_s h"
	
	return
	end subroutine
