	subroutine eval_graphite(iaq, ibq)

	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, ibq, ib4, ib216, ib64to216, ib512, ib64to512, ib1000, ib64to1000
	integer*4 i, j, ictk, mi, ni
	real*8 sgcw6, epcw, xi, yi, zi, xij, yij, zij, rijsq, rij, rsi, r6, rpl, flj
	
	sgcw6 = (3.19D0 / sigma)**6
	epcw = 392.0D0 / epsilon
	
	if(ibq.eq.0) then
	
		do ib4 = 1, 16
		
			ib216 = ib64to216(ib4)
			ib512 = ib64to512(ib4)
			ib1000 = ib64to1000(ib4)
			
			ictk = 3 * ictw1(ib216)
			do i = 1, ictk - 2, 3
			
				xi = rxw1(i, ib216)
				yi = ryw1(i, ib216)
				zi = rzw1(i, ib216)
				
				do j = 1, ictgraphite(ib4)
				
					xij = xi - xgr(j, ib4)
					yij = yi - ygr(j, ib4)
					zij = zi - zgr(j, ib4)
					
					rijsq = xij**2 + yij**2 + zij**2
					rij = sqrt(rijsq)
					
					if(rij.le.rc) then
					
						rsi = 1.0D0 / rijsq
						r6 = rsi**3
						rpl = 48.0D0 * epcw * sgcw6 * r6 * (sgcw6 * r6 - 0.5D0)
						flj = rpl * rsi
						
						fxw1(i, ib216) = fxw1(i, ib216) + xij * flj
						fyw1(i, ib216) = fyw1(i, ib216) + yij * flj
						fzw1(i, ib216) = fzw1(i, ib216) + zij * flj
						
					endif
					
				enddo
				
			enddo
			
			ictk = ictW2Central(ib4)
			do i = 1, ictk
			
				mi = W2Central(i, 3, ib4)
				ni = W2Central(i, 4, ib4)
				
				if(ni.eq.1) then
				
					xi = rxw2(mi, ib512)
					yi = ryw2(mi, ib512)
					zi = rzw2(mi, ib512)
					
					do j = 1, ictgraphite(ib4)
					
						xij = xi - xgr(j, ib4)
						yij = yi - ygr(j, ib4)
						zij = zi - zgr(j, ib4)
						
						rijsq = xij**2 + yij**2 + zij**2
						rij = sqrt(rijsq)
						
						if(rij.le.rc) then
						
							rsi = 1.0D0 / rijsq
							r6 = rsi**3
							rpl = 48.0D0 * epcw * sgcw6 * r6 * (sgcw6 * r6 - 0.5D0)
							flj = rpl * rsi
							
							fxw2(mi, ib512) = fxw2(mi, ib512) + xij * flj
							fyw2(mi, ib512) = fyw2(mi, ib512) + yij * flj
							fzw2(mi, ib512) = fzw2(mi, ib512) + zij * flj
							
						endif
						
					enddo
					
				endif
				
			enddo
			
			ictk = ictSCentral(ib4)
			do i = 1, ictk
			
				mi = SCentral(i, 3, ib4)
				ni = SCentral(i, 4, ib4)
				
				xi = rxs(mi, ib1000)
				yi = rys(mi, ib1000)
				zi = rzs(mi, ib1000)
				
				do j = 1, ictgraphite(ib4)
				
					xij = xi - xgr(j, ib4)
					yij = yi - ygr(j, ib4)
					zij = zi - zgr(j, ib4)
					
					rijsq = xij**2 + yij**2 + zij**2
					rij = sqrt(rijsq)
					
					if(rij.le.rc) then
					
						rsi = 1.0D0 / rijsq
						r6 = rsi**3
						rpl = 48.0D0 * epcs(ni) * sgcs6(ni) * r6 * (sgcs6(ni) * r6 - 0.5D0)
						flj = rpl * rsi
						
						fxs(mi, ib1000) = fxs(mi, ib1000) + xij * flj
						fys(mi, ib1000) = fys(mi, ib1000) + yij * flj
						fzs(mi, ib1000) = fzs(mi, ib1000) + zij * flj
						!write(*,*) iaq, ibq, epcs(ni), sgcs6(ni), ni
					endif
					
				enddo
				
			enddo
		
		enddo
	
	endif
	
	return
	end subroutine
