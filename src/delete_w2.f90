	subroutine delete_w2(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, iaq1
	integer*4 i, j, k, l, g, h
	integer*4 ib4, ib216, ib64to216, ib512, ib512m, ib64to512, ibox2, ibox64, ibox64mod
	integer*4 ibO, ibG, ibH
	integer*4 ib64O, ib64G, ib64H
	integer*4 ib64Om, ib64Gm, ib64Hm
	integer*4 ictk, istart, jstart
	integer*4 ib, ibm, ib64, ib64m, ib216m, iwlist
	
	real*8 xO, yO, zO
	real*8 xG, yG, zG
	real*8 xH, yH, zH
	real*8 rxL2, ryL2, rzL2
	real*8 rx216m, ry216m, rz216m

	real*8 rwatx(kdeletew23), rwaty(kdeletew23), rwatz(kdeletew23)
	real*8 vwatx(kdeletew23), vwaty(kdeletew23), vwatz(kdeletew23)

	iaq1 = iaq + 1
	
	jdeletedw2 = 0
	
	g = 0
	h = 0
	
	do ib4 = 1, 64
	
		k = 0
		
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
			
			if(ibO.eq.iaq1) then
			
				if(ibO.eq.ibG.and.ibO.eq.ibH) then
				
					ib64O = ibox64(xO, yO, zO)
					ib64G = ibox64(xG, yG, zG)
					ib64H = ibox64(xH, yH, zH)
					
					if(ib64O.eq.ib64G.and.ib64O.eq.ib64H) then
					!write(*,*) "w2->w1"
						ib216m = ib64to216(ib64O)
						istart = ictw1(ib216m) * 3
					
						do l = 1, 3
						
							rxw1(istart + l, ib216m) = rxw2(jstart + l, ib512)
							ryw1(istart + l, ib216m) = ryw2(jstart + l, ib512)
							rzw1(istart + l, ib216m) = rzw2(jstart + l, ib512)
							
							vxw1(istart + l, ib64O) = vxw2(jstart + l, ib4)
							vyw1(istart + l, ib64O) = vyw2(jstart + l, ib4)
							vzw1(istart + l, ib64O) = vzw2(jstart + l, ib4)
						
						enddo
						
						ictw1(ib216m) = ictw1(ib216m) + 1
						
					else
					
						if(ib64O.eq.ib4) then
						
							do l = 1, 3
							
								k = k + 1
								rwatx(k) = rxw2(jstart + l, ib512)
								rwaty(k) = ryw2(jstart + l, ib512)
								rwatz(k) = rzw2(jstart + l, ib512)
								
								vwatx(k) = vxw2(jstart + l, ib4)
								vwaty(k) = vyw2(jstart + l, ib4)
								vwatz(k) = vzw2(jstart + l, ib4)
								
							enddo
						
						else
						!write(*,*) "w2->w2"
							ib512m = ib64to512(ib64O)
							istart = ictw2(ib512m) * 3
							
							do l = 1, 3
							
								rxw2(istart + l, ib512m) = rxw2(jstart + l, ib512)
								ryw2(istart + l, ib512m) = ryw2(jstart + l, ib512)
								rzw2(istart + l, ib512m) = rzw2(jstart + l, ib512)
								
								vxw2(istart + l, ib64O) = vxw2(jstart + l, ib4)
								vyw2(istart + l, ib64O) = vyw2(jstart + l, ib4)
								vzw2(istart + l, ib64O) = vzw2(jstart + l, ib4)
								
							enddo
							
							ictw2(ib512m) = ictw2(ib512m) + 1
							
						endif
					
					endif
				
				! ibO equals iaq1 but ibO not equal to both ibG and ibH
				else
				
					ib64O = ibox64(xO, yO, zO)
				
					if(ib64O.eq.ib4) then
						
						do l = 1, 3
						
							k = k + 1
							rwatx(k) = rxw2(jstart + l, ib512)
							rwaty(k) = ryw2(jstart + l, ib512)
							rwatz(k) = rzw2(jstart + l, ib512)
							
							vwatx(k) = vxw2(jstart + l, ib4)
							vwaty(k) = vyw2(jstart + l, ib4)
							vwatz(k) = vzw2(jstart + l, ib4)
							
						enddo
					
					else
					!write(*,*) "w2->w2"
						ib512m = ib64to512(ib64O)
						istart = ictw2(ib512m) * 3
						
						do l = 1, 3
						
							rxw2(istart + l, ib512m) = rxw2(jstart + l, ib512)
							ryw2(istart + l, ib512m) = ryw2(jstart + l, ib512)
							rzw2(istart + l, ib512m) = rzw2(jstart + l, ib512)
							
							vxw2(istart + l, ib64O) = vxw2(jstart + l, ib4)
							vyw2(istart + l, ib64O) = vyw2(jstart + l, ib4)
							vzw2(istart + l, ib64O) = vzw2(jstart + l, ib4)
							
						enddo
						
						ictw2(ib512m) = ictw2(ib512m) + 1
					
					endif
					
				endif
			
			! key atom of molecule not in iaq
			else
			
				jdeletedw2 = jdeletedw2 + 1
			
				rxL2 = rcen2x(ibO)
				ryL2 = rcen2y(ibO)
				rzL2 = rcen2z(ibO)
			
				ib64Om = ibox64mod(xO, yO, zO, rxL2, ryL2, rzL2)
			
				if(ibO.eq.ibG.and.ibO.eq.ibH) then
				
					ib64Gm = ibox64mod(xG, yG, zG, rxL2, ryL2, rzL2)
					ib64Hm = ibox64mod(xH, yH, zH, rxL2, ryL2, rzL2)
					
					if(ib64Om.eq.ib64Gm.and.ib64Om.eq.ib64Hm) then
						iwlist = 1
					else
						iwlist = 2
					endif
					
				else
				
					iwlist = 2
					
				endif
				
				iw2ewnsud(g + 1) = ibO - 1
				iw2ewnsud(g + 2) = ib64Om
				iw2ewnsud(g + 3) = iwlist
				g = g + 3
				
				do l = 1, 3
				
					rw2ewnsud(h + 1) = rxw2(jstart + l, ib512)
					rw2ewnsud(h + 2) = ryw2(jstart + l, ib512)
					rw2ewnsud(h + 3) = rzw2(jstart + l, ib512)
					
					rw2ewnsud(h + 4) = vxw2(jstart + l, ib4)
					rw2ewnsud(h + 5) = vyw2(jstart + l, ib4)
					rw2ewnsud(h + 6) = vzw2(jstart + l, ib4)
					h = h + 6
					
				enddo
			
			endif
			
		enddo
		
		! assign new list
		ictw2(ib512) = k / 3
		do i = 1, k
		
			rxw2(i, ib512) = rwatx(i)
			ryw2(i, ib512) = rwaty(i)
			rzw2(i, ib512) = rwatz(i)
			
			vxw2(i, ib4) = vwatx(i)
			vyw2(i, ib4) = vwaty(i)
			vzw2(i, ib4) = vwatz(i)
			
		enddo

		!if(ictk.ne.k/3) write(*,'(i5,i5,i5,i5,i5)') iaq, ib4, ictk, k / 3, 2

		if(k.gt.kdeletew23) write(*,*) "delete_w2 k"
		if(ictw2(ib512).gt.nw) write(*,*) "delete_w2 ictw2"
	
	enddo

	!if(jdeletedw2.ne.0) write(*,*) "w2", iaq, jdeletedw2
	if(h.gt.kdeletew236) write(*,*) "delete_w2 h"
	
	return
	end subroutine
