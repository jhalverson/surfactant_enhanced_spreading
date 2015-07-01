	subroutine eval_wall()

	implicit none
	include "mybox.cmns"
	
	integer*4 i, ib4, ib216, ib64to216, ib512, ib64to512, ib1000, ib64to1000, ictk

	do ib4 = 1, 64

		ib216 = ib64to216(ib4)
		ib512 = ib64to512(ib4)
		ib1000 = ib64to1000(ib4)

		ictk = 3 * ictw1(ib216)
		do i = 1, ictk

			fxw1(i, ib216) = fxw1(i, ib216) + 48.0D0 / (rxw1(i, ib216) - rhsx)**13
			fxw1(i, ib216) = fxw1(i, ib216) + 48.0D0 / (rhsx + rxw1(i, ib216))**13
			
			fyw1(i, ib216) = fyw1(i, ib216) + 48.0D0 / (ryw1(i, ib216) - rhsy)**13
			fyw1(i, ib216) = fyw1(i, ib216) + 48.0D0 / (rhsy + ryw1(i, ib216))**13
			
			fzw1(i, ib216) = fzw1(i, ib216) + 48.0D0 / (rzw1(i, ib216) - rhsz)**13
			fzw1(i, ib216) = fzw1(i, ib216) + 48.0D0 / (rhsz + rzw1(i, ib216))**13
			
		enddo
		
		ictk = 3 * ictw2(ib512)
		do i = 1, ictk

			fxw2(i, ib512) = fxw2(i, ib512) + 48.0D0 / (rxw2(i, ib512) - rhsx)**13
			fxw2(i, ib512) = fxw2(i, ib512) + 48.0D0 / (rhsx + rxw2(i, ib512))**13
			
			fyw2(i, ib512) = fyw2(i, ib512) + 48.0D0 / (ryw2(i, ib512) - rhsy)**13
			fyw2(i, ib512) = fyw2(i, ib512) + 48.0D0 / (rhsy + ryw2(i, ib512))**13
			
			fzw2(i, ib512) = fzw2(i, ib512) + 48.0D0 / (rzw2(i, ib512) - rhsz)**13
			fzw2(i, ib512) = fzw2(i, ib512) + 48.0D0 / (rhsz + rzw2(i, ib512))**13
			
		enddo
		
		ictk = 29 * icts(ib1000)
		do i = 1, ictk

			fxs(i, ib1000) = fxs(i, ib1000) + 48.0D0 / (rxs(i, ib1000) - rhsx)**13
			fxs(i, ib1000) = fxs(i, ib1000) + 48.0D0 / (rhsx + rxs(i, ib1000))**13
			
			fys(i, ib1000) = fys(i, ib1000) + 48.0D0 / (rys(i, ib1000) - rhsy)**13
			fys(i, ib1000) = fys(i, ib1000) + 48.0D0 / (rhsy + rys(i, ib1000))**13
			
			fzs(i, ib1000) = fzs(i, ib1000) + 48.0D0 / (rzs(i, ib1000) - rhsz)**13
			fzs(i, ib1000) = fzs(i, ib1000) + 48.0D0 / (rhsz + rzs(i, ib1000))**13
			
		enddo
	
	enddo
	
	return
	end subroutine
