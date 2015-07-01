	subroutine eval_surfactant_internal(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 i, j, k, m, n, iaq, ib4, ictk, ib1000, ib64to1000, isf, jsf
	integer*4 jstart, jend, kstart, kend, jtype, ktype
	
	real*8 xj, yj, zj, x, y, z
	real*8 rsq, rctr, fel, flj, ft, rsi, r6, rpl
	
	do ib4 = 1, 64
	
		ib1000 = ib64to1000(ib4)
		
		ictk = icts(ib1000)
		do i = 1, ictk
		
			jstart = (i - 1) * 29 + 1
			jend = (i - 1) * 29 + 15
			
			do j = jstart, jend
			
				jtype = j - (i - 1) * 29
				
				xj = rxs(j, ib1000)
				yj = rys(j, ib1000)
				zj = rzs(j, ib1000)
				
				kstart = j + 4
				kend = (i - 1) * 29 + 29
				
				do k = kstart, kend
				
					ktype = k - (i - 1) * 29
					
					x = xj - rxs(k, ib1000)
					y = yj - rys(k, ib1000)
					z = zj - rzs(k, ib1000)
					
					rsq = x * x + y * y + z * z
					rctr = sqrt(rsq)
					
					fel = rqs(jtype) * rqs(ktype) * xpi / rctr**3
					!if(j.eq.1) write(*,*) jtype, ktype, fel
					flj = 0.0D0
					!if(rctr.le.rc) then
					
						rsi = 1.0D0 / rsq
						r6 = rsi**3
						rpl = 48.0D0 * epss(jtype, ktype) * sgss6(jtype, ktype) * r6 * (sgss6(jtype, ktype) * r6 - 0.5D0)
						flj = rpl * rsi
						
					!endif
					!if(iaq.eq.14.and.ib1000.eq.556) write(*,*) j, k, flj, fxs(1, ib1000)
					ft = fel + flj
					
					fxs(j, ib1000) = fxs(j, ib1000) + x * ft
					fys(j, ib1000) = fys(j, ib1000) + y * ft
					fzs(j, ib1000) = fzs(j, ib1000) + z * ft
					
					fxs(k, ib1000) = fxs(k, ib1000) - x * ft
					fys(k, ib1000) = fys(k, ib1000) - y * ft
					fzs(k, ib1000) = fzs(k, ib1000) - z * ft
					
				enddo
				
			enddo
			
			do k = 1, 42
			
				m = iself(k)
				n = jself(k)
				
				isf = (i - 1) * 29 + m
				jsf = (i - 1) * 29 + n
				
				x = rxs(isf, ib1000) - rxs(jsf, ib1000)
				y = rys(isf, ib1000) - rys(jsf, ib1000)
				z = rzs(isf, ib1000) - rzs(jsf, ib1000)
				
				rsq = x * x + y * y + z * z
				rctr = sqrt(rsq)
				
				fel = rqs(m) * rqs(n) * xpi / rctr**3
				!write(*,*) m, n, fel
				flj = 0.0D0
				!if(rctr.le.rc) then
				
					rsi = 1.0D0 / rsq
					r6 = rsi**3
					rpl = 48.0D0 * epss(m, n) * sgss6(m, n) * r6 * (sgss6(m, n) * r6 - 0.5D0)
					flj = rpl * rsi
				
				!endif
				!if(iaq.eq.14.and.ib1000.eq.556) write(*,*) isf, jsf, flj, fxs(1, ib1000)
				
				ft = fel + flj
				
				fxs(isf, ib1000) = fxs(isf, ib1000) + x * ft
				fys(isf, ib1000) = fys(isf, ib1000) + y * ft
				fzs(isf, ib1000) = fzs(isf, ib1000) + z * ft
				
				fxs(jsf, ib1000) = fxs(jsf, ib1000) - x * ft
				fys(jsf, ib1000) = fys(jsf, ib1000) - y * ft
				fzs(jsf, ib1000) = fzs(jsf, ib1000) - z * ft
				
			enddo
			
		enddo
		
	enddo
	
	return
	end subroutine
