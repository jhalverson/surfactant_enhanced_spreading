	subroutine delete_w1(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, iaq1
	integer*4 i, j, k, l, g, h
	integer*4 ib4, ib216, ib64to216, ib512m, ib64to512, ibox2, ibox64, ibox64mod
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
	
	real*8 rwatx(kdeletew13), rwaty(kdeletew13), rwatz(kdeletew13)
	real*8 vwatx(kdeletew13), vwaty(kdeletew13), vwatz(kdeletew13)

	iaq1 = iaq + 1
	
	jdeletedw1 = 0
	
	g = 0
	h = 0
	
	do ib4 = 1, 64
	
		k = 0
		
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
			
			if(ibO.eq.iaq1) then
			
				if(ibO.eq.ibG.and.ibO.eq.ibH) then
				
					ib64O = ibox64(xO, yO, zO)
					ib64G = ibox64(xG, yG, zG)
					ib64H = ibox64(xH, yH, zH)
					
					if(ib64O.eq.ib64G.and.ib64O.eq.ib64H) then
					
						if(ib64O.eq.ib4) then
						
							do l = 1, 3
							
								k = k + 1
								rwatx(k) = rxw1(jstart + l, ib216)
								rwaty(k) = ryw1(jstart + l, ib216)
								rwatz(k) = rzw1(jstart + l, ib216)
								
								vwatx(k) = vxw1(jstart + l, ib4)
								vwaty(k) = vyw1(jstart + l, ib4)
								vwatz(k) = vzw1(jstart + l, ib4)
								
							enddo
						
						else
							!write(*,*) "unlikely here w1; iaq =", iaq
							ib216m = ib64to216(ib64O)
							istart = ictw1(ib216m) * 3
							
							do l = 1, 3
							
								rxw1(istart + l, ib216m) = rxw1(jstart + l, ib216)
								ryw1(istart + l, ib216m) = ryw1(jstart + l, ib216)
								rzw1(istart + l, ib216m) = rzw1(jstart + l, ib216)
								
								vxw1(istart + l, ib64O) = vxw1(jstart + l, ib4)
								vyw1(istart + l, ib64O) = vyw1(jstart + l, ib4)
								vzw1(istart + l, ib64O) = vzw1(jstart + l, ib4)
							
							enddo
							
							ictw1(ib216m) = ictw1(ib216m) + 1
						
						endif
					
					! ibO equals ibG and ibO equals ibH but ib64O not equal to both ib64G and ib64H
					else
						!write(*,*) "w1->w2", ib4, ib64O, ib64G, ib64H
						ib512m = ib64to512(ib64O)
						istart = ictw2(ib512m) * 3
						
						do l = 1, 3
						
							rxw2(istart + l, ib512m) = rxw1(jstart + l, ib216)
							ryw2(istart + l, ib512m) = ryw1(jstart + l, ib216)
							rzw2(istart + l, ib512m) = rzw1(jstart + l, ib216)
							
							vxw2(istart + l, ib64O) = vxw1(jstart + l, ib4)
							vyw2(istart + l, ib64O) = vyw1(jstart + l, ib4)
							vzw2(istart + l, ib64O) = vzw1(jstart + l, ib4)
							
						enddo
						
						ictw2(ib512m) = ictw2(ib512m) + 1
					
					endif
				
				! ibO equals iaq1 but ibO not equal to both ibG and ibH
				else
					!write(*,*) "at least 1 H out of iaq"
					ib64O = ibox64(xO, yO, zO)
					ib512m = ib64to512(ib64O)
					istart = ictw2(ib512m) * 3
					
					do l = 1, 3
					
						rxw2(istart + l, ib512m) = rxw1(jstart + l, ib216)
						ryw2(istart + l, ib512m) = ryw1(jstart + l, ib216)
						rzw2(istart + l, ib512m) = rzw1(jstart + l, ib216)
						
						vxw2(istart + l, ib64O) = vxw1(jstart + l, ib4)
						vyw2(istart + l, ib64O) = vyw1(jstart + l, ib4)
						vzw2(istart + l, ib64O) = vzw1(jstart + l, ib4)
						
					enddo
					
					ictw2(ib512m) = ictw2(ib512m) + 1
					
				endif
			
			! key atom of molecule not in iaq
			else
			
				jdeletedw1 = jdeletedw1 + 1
				
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
				
				iw1ewnsud(g + 1) = ibO - 1
				iw1ewnsud(g + 2) = ib64Om
				iw1ewnsud(g + 3) = iwlist
				g = g + 3
				
				do l = 1, 3
				
					rw1ewnsud(h + 1) = rxw1(jstart + l, ib216)
					rw1ewnsud(h + 2) = ryw1(jstart + l, ib216)
					rw1ewnsud(h + 3) = rzw1(jstart + l, ib216)
					
					rw1ewnsud(h + 4) = vxw1(jstart + l, ib4)
					rw1ewnsud(h + 5) = vyw1(jstart + l, ib4)
					rw1ewnsud(h + 6) = vzw1(jstart + l, ib4)
					h = h + 6
					
				enddo
			
			endif
			
		enddo
		
		! assign new list
		ictw1(ib216) = k / 3
		do i = 1, k
		
			rxw1(i, ib216) = rwatx(i)
			ryw1(i, ib216) = rwaty(i)
			rzw1(i, ib216) = rwatz(i)
			
			vxw1(i, ib4) = vwatx(i)
			vyw1(i, ib4) = vwaty(i)
			vzw1(i, ib4) = vwatz(i)
		
		enddo

		!if(ictk.ne.k/3.and.iaq.eq.14) write(*,'(i5,i5,i5,i5,i5)') iaq, ib4, ictk, k / 3, 1

		if(k.gt.kdeletew13) write(*,*) "delete_w1 k"
		if(ictw1(ib216).gt.nw) write(*,*) "delete_w1 ictw1"
	
	enddo
	
	!if(jdeletedw1.ne.0) write(*,*) "w1", iaq, jdeletedw1
	if(h.gt.kdeletew136) write(*,*) "delete_w1 h"
	
	return
	end subroutine
