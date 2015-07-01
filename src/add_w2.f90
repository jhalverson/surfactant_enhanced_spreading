	subroutine add_w2(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq
	integer*4 ib1, ib64, ib216, ib64to216, ib512, ib64to512
	integer*4 l, j, jstart, h, iwlist
	
	do j = 1, jdeletedw2_total
	
		ib1 = iw2ewnsud(3 * j - 2)
		
		if(ib1.eq.iaq) then
		
			ib64 = iw2ewnsud(3 * j - 1)
			iwlist = iw2ewnsud(3 * j)
			
			if(iwlist.eq.1) then
			
				ib216 = ib64to216(ib64)
				jstart = 3 * ictw1(ib216)
				h = (j - 1) * 18
				
				do l = 1, 3
				
					rxw1(jstart + l, ib216) = rw2ewnsud(h + 1)
					ryw1(jstart + l, ib216) = rw2ewnsud(h + 2)
					rzw1(jstart + l, ib216) = rw2ewnsud(h + 3)
				
					vxw1(jstart + l, ib64) = rw2ewnsud(h + 4)
					vyw1(jstart + l, ib64) = rw2ewnsud(h + 5)
					vzw1(jstart + l, ib64) = rw2ewnsud(h + 6)
					
					h = h + 6
				
				enddo
				
				ictw1(ib216) = ictw1(ib216) + 1
				if(ictw1(ib216).gt.nw) write(*,*) "add_w2 nw w1"
				!write(*,*) "add to w1 (add_w2)"
			
			else
			
				ib512 = ib64to512(ib64)
				jstart = 3 * ictw2(ib512)
				h = (j - 1) * 18
				
				do l = 1, 3
			
					rxw2(jstart + l, ib512) = rw2ewnsud(h + 1)
					ryw2(jstart + l, ib512) = rw2ewnsud(h + 2)
					rzw2(jstart + l, ib512) = rw2ewnsud(h + 3)
					
					vxw2(jstart + l, ib64) = rw2ewnsud(h + 4)
					vyw2(jstart + l, ib64) = rw2ewnsud(h + 5)
					vzw2(jstart + l, ib64) = rw2ewnsud(h + 6)
					
					h = h + 6
	
				enddo
				
				ictw2(ib512) = ictw2(ib512) + 1
				!write(*,*) ictw2(ib512)
				if(ictw2(ib512).gt.nw) write(*,*) "add_w2 nw w2"
				!write(*,*) "add to w2 (add_w2)"
			
			endif
			
		endif
	
	enddo
	
	return
	end subroutine
