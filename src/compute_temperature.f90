	subroutine compute_temperature(iaq)
	
	implicit none
	include "mpif.h"
	include "mybox.cmns"
	
	integer*4 iaq, ib4, ib216, ib64to216, ib512, ib64to512, ib1000, ib64to1000
	integer*4 i, ictk, ierr
	
	real*8 sumvelO, sumvelH, sumvelCH2, sumvelSI, sumvelCH3
	real*8 v2sum, v2total, temp
	
	sumvelO = 0.0D0
	sumvelH = 0.0D0
	sumvelCH2 = 0.0D0
	sumvelSI = 0.0D0
	sumvelCH3 = 0.0D0
	
	do ib4 = 1, 64
	
		ib216 = ib64to216(ib4)
		ib512 = ib64to512(ib4)
		ib1000 = ib64to1000(ib4)
		
		ictk = 3 * ictw1(ib216)
		do i = 1, ictk - 2, 3
		
			sumvelO = sumvelO + vxw1(i, ib4)**2 + vyw1(i, ib4)**2 + vzw1(i, ib4)**2
			
		enddo
		
		do i = 2, ictk - 1, 3
		
			sumvelH = sumvelH + vxw1(i, ib4)**2 + vyw1(i, ib4)**2 + vzw1(i, ib4)**2
			
		enddo
		
		do i = 3, ictk - 0, 3
		
			sumvelH = sumvelH + vxw1(i, ib4)**2 + vyw1(i, ib4)**2 + vzw1(i, ib4)**2
			
		enddo
		
		ictk = 3 * ictw2(ib512)
		do i = 1, ictk - 2, 3
		
			sumvelO = sumvelO + vxw2(i, ib4)**2 + vyw2(i, ib4)**2 + vzw2(i, ib4)**2
			
		enddo
		
		do i = 2, ictk - 1, 3
		
			sumvelH = sumvelH + vxw2(i, ib4)**2 + vyw2(i, ib4)**2 + vzw2(i, ib4)**2
			
		enddo
		
		do i = 3, ictk - 0, 3
		
			sumvelH = sumvelH + vxw2(i, ib4)**2 + vyw2(i, ib4)**2 + vzw2(i, ib4)**2
			
		enddo
		
		ictk = 29 * icts(ib1000)
		do i = 1, ictk - 28, 29
		
			sumvelH = sumvelH + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 2, ictk - 27, 29
		
			sumvelO = sumvelO + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 3, ictk - 26, 29
		
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 4, ictk - 25, 29
		
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 5, ictk - 24, 29
		
			sumvelO = sumvelO + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 6, ictk - 23, 29
		
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 7, ictk - 22, 29
		
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 8, ictk - 21, 29
		
			sumvelO = sumvelO + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 9, ictk - 20, 29
		
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 10, ictk - 19, 29
		
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 11, ictk - 18, 29
		
			sumvelO = sumvelO + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 12, ictk - 17, 29
		
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 13, ictk - 16, 29
		
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 14, ictk - 15, 29
		
			sumvelO = sumvelO + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 15, ictk - 14, 29
		
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 16, ictk - 13, 29
		
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 17, ictk - 12, 29
		
			sumvelCH2 = sumvelCH2 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 18, ictk - 11, 29
		
			sumvelSI = sumvelSI + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 19, ictk - 10, 29
		
			sumvelCH3 = sumvelCH3 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 20, ictk - 9, 29
		
			sumvelO = sumvelO + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 21, ictk - 8, 29
		
			sumvelSI = sumvelSI + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 22, ictk - 7, 29
		
			sumvelCH3 = sumvelCH3 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 23, ictk - 6, 29
		
			sumvelCH3 = sumvelCH3 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 24, ictk - 5, 29
		
			sumvelCH3 = sumvelCH3 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 25, ictk - 4, 29
		
			sumvelO = sumvelO + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 26, ictk - 3, 29
		
			sumvelSI = sumvelSI + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 27, ictk - 2, 29
		
			sumvelCH3 = sumvelCH3 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 28, ictk - 1, 29
		
			sumvelCH3 = sumvelCH3 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
			
		enddo
		
		do i = 29, ictk - 0, 29
		
			sumvelCH3 = sumvelCH3 + vxs(i, ib4)**2 + vys(i, ib4)**2 + vzs(i, ib4)**2
		
		enddo
	
	enddo
	
	v2sum = emw(1) * sumvelO + emw(2) * sumvelH + ems(3) * sumvelCH2 + ems(18) * sumvelSI + ems(19) * sumvelCH3
	call mpi_allreduce(v2sum, v2total, 1, mpi_real8, mpi_sum, mpi_comm_world, ierr)
	temp = v2total / dtsq / (3 * (3 * iwater + 29 * isurfactant) - (3 * iwater + 28 * isurfactant))
	
	if(iaq.eq.0) write(*,*) "T =", temp
	
	return
	end subroutine
