	subroutine add_s(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq
	integer*4 ib1, ib64, ib1000, ib64to1000
	integer*4 l, j, jstart, h
	
	do j = 1, jdeleteds_total
	
		ib1 = isewnsud(2 * j - 1)
		
		if(ib1.eq.iaq) then
		
			ib64 = isewnsud(2 * j)
			ib1000 = ib64to1000(ib64)
			
			jstart = 29 * icts(ib1000)
			h = (j - 1) * 174
			
			do l = 1, 29
			
				rxs(jstart + l, ib1000) = rsewnsud(h + 1)
				rys(jstart + l, ib1000) = rsewnsud(h + 2)
				rzs(jstart + l, ib1000) = rsewnsud(h + 3)

				vxs(jstart + l, ib64) = rsewnsud(h + 4)
				vys(jstart + l, ib64) = rsewnsud(h + 5)
				vzs(jstart + l, ib64) = rsewnsud(h + 6)

				h = h + 6

			enddo

			icts(ib1000) = icts(ib1000) + 1
			if(icts(ib1000).gt.ns) write(*,*) "add_s ns"

		endif

	enddo

	return
	end subroutine
