	subroutine read_graphite(iaq, ibq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, ibq, i, ib4
	real*8 xtmp, ytmp, ztmp, rcut, rsq, rctr4x, rctr4y
	
	rcut = (rc + sqrt(2.0D0) * rs / 2.0D0)**2
	ictgraphite = 0
	
	if(ibq.eq.0) then
	
		do ib4 = 1, 16
		
			rctr4x = rcen4x(ib4)
			rctr4y = rcen4y(ib4)
			
			open(10, file = "graphite.dat")
			do i = 1, 22956
			
				read(10,'(f10.4,f11.4,f11.4)') xtmp, ytmp, ztmp
				rsq = (xtmp - rctr4x)**2 + (ytmp - rctr4y)**2
				
				if(rsq.le.rcut) then
				
					ictgraphite(ib4) = ictgraphite(ib4) + 1
					if(ictgraphite(ib4).ge.kgraphite) write(*,*) "graphite bound exceeded."
					xgr(ictgraphite(ib4), ib4) = xtmp
					ygr(ictgraphite(ib4), ib4) = ytmp
					zgr(ictgraphite(ib4), ib4) = ztmp
					
				endif
				
			enddo
			close(10)
			
		enddo
	
	endif
	
	return
	end subroutine
