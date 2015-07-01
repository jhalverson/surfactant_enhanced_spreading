	subroutine move_w1_b(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, ib216, ib64to216
	integer*4 nb, i, j, ictmolecule
	integer*4 iatom, ib4, it, ib, jb, maxit
	real*8 tol, tol2
	real*8 vxab, vyab, vzab, rvab
	real*8 rxab, ryab, rzab, rma, rmb, gab
	real*8 dx, dy, dz
	real*8 rxi(3), ryi(3), rzi(3)
	real*8 vxi(3), vyi(3), vzi(3)
	real*8 escale
	logical ldone
	
	nb = 3
	
	tol = 0.0001D0
	maxit = 100
	
	escale = 2.0D0 / (2.0D0 + eta * dt)
	
	do ib4 = 1, 64
	
		ib216 = ib64to216(ib4)
		
		ictmolecule = ictw1(ib216)
		do j = 1, ictmolecule
		
			do i = 1, 3
			
			iatom = 3 * (j - 1) + i
			
			rxi(i) = rxw1(iatom, ib216)
			ryi(i) = ryw1(iatom, ib216)
			rzi(i) = rzw1(iatom, ib216)
			
			axw1(iatom, ib4) = fxw1(iatom, ib216) * dtsq2 / emw(i)
			ayw1(iatom, ib4) = fyw1(iatom, ib216) * dtsq2 / emw(i)
			azw1(iatom, ib4) = fzw1(iatom, ib216) * dtsq2 / emw(i)
			
			vxi(i) = escale * (vxw1(iatom, ib4) + axw1(iatom, ib4))
			vyi(i) = escale * (vyw1(iatom, ib4) + ayw1(iatom, ib4))
			vzi(i) = escale * (vzw1(iatom, ib4) + azw1(iatom, ib4))
			
		enddo
		   
		it = 0
		ldone = .false.
		   
1009	   if((.not.ldone).and.(it.le.maxit)) then
			  
			  ldone = .true.
			  
			  do i = 1, nb
				 
				 ib = i + 1
				 if(ib.gt.3) ib = 1
				 
				 vxab = vxi(i) - vxi(ib)
				 vyab = vyi(i) - vyi(ib)
				 vzab = vzi(i) - vzi(ib)
				 
				 rxab = rxi(i) - rxi(ib)
				 ryab = ryi(i) - ryi(ib)
				 rzab = rzi(i) - rzi(ib)
				 
				 rvab = rxab * vxab + ryab * vyab + rzab * vzab
				 rvab = rvab / dt
				 rma = 1.0D0 / emw(i)
				 rmb = 1.0D0 / emw(ib)
				 gab = -rvab / ((rma + rmb) * dsqw(i))
				 
				 if(abs(gab).gt.tol) then
					
					dx = rxab * gab
					dy = ryab * gab
					dz = rzab * gab
					
					vxi(i) = vxi(i) + rma * dx * dt
					vyi(i) = vyi(i) + rma * dy * dt
					vzi(i) = vzi(i) + rma * dz * dt
					
					vxi(ib) = vxi(ib) - rmb * dx * dt
					vyi(ib) = vyi(ib) - rmb * dy * dt
					vzi(ib) = vzi(ib) - rmb * dz * dt
					
					ldone = .false.
					
				 endif
				 
			  enddo
			  
			  it = it + 1
			  goto 1009
			  
		   endif
		   
		   if(.not.ldone) then
			  
			  write(*,*) "iaq", iaq, j, "Too many iterations in move_water_b."
			  open(10, file = "ERROR.txt")
			  write(10,*) j, "Too many iterations in move_water_b."
			  close(10)
			  !stop
			  
		   endif
		   
		   do i = 1, 3
			  
			  iatom = 3 * (j - 1) + i
			  
			  vxw1(iatom, ib4) = vxi(i)
			  vyw1(iatom, ib4) = vyi(i)
			  vzw1(iatom, ib4) = vzi(i)
			  
		   enddo
		   
		enddo
	
	enddo
	
	return
	end subroutine
