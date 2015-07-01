	subroutine move_w1_a(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, ib216, ib64to216
	integer*4 nb, i, j, ictmolecule
	integer*4 iatom, ib4, it, ib, jb, maxit
	real*8 tol, tol2
	real*8 pxab, pyab, pzab, pabsq, rabsq, diffsq
	real*8 rxab, ryab, rzab, rpab, rma, rmb, gab
	real*8 dx, dy, dz
	real*8 rxi(3), ryi(3), rzi(3)
	real*8 pxi(3), pyi(3), pzi(3)
	real*8 vxi(3), vyi(3), vzi(3)
	real*8 escale
	logical ldone
	
	tol = 0.0001D0
	maxit = 100
	
	nb = 3
	tol2 = 2.0D0 * tol
	
	escale = 1.0D0 - eta * dt / 2.0D0
	
	call nh_water1(w1v1)
	
	do ib4 = 1, 64
	
		ib216 = ib64to216(ib4)
		
		ictmolecule = ictw1(ib216)
		do j = 1, ictmolecule
		
		do i = 1, 3
			
			iatom = (j - 1) * 3 + i
			
			rxi(i) = rxw1(iatom, ib216)
			ryi(i) = ryw1(iatom, ib216)
			rzi(i) = rzw1(iatom, ib216)
			
			pxi(i) = rxw1(iatom, ib216) + vxw1(iatom, ib4) * escale + axw1(iatom, ib4)
			pyi(i) = ryw1(iatom, ib216) + vyw1(iatom, ib4) * escale + ayw1(iatom, ib4)
			pzi(i) = rzw1(iatom, ib216) + vzw1(iatom, ib4) * escale + azw1(iatom, ib4)
			
			vxi(i) = vxw1(iatom, ib4) * escale + axw1(iatom, ib4)
			vyi(i) = vyw1(iatom, ib4) * escale + ayw1(iatom, ib4)
			vzi(i) = vzw1(iatom, ib4) * escale + azw1(iatom, ib4)
			  
		enddo
		   
		it = 0
		ldone = .false.
		   
	 1000  if((.not.ldone).and.(it.le.maxit)) then
		   ldone = .true.
			  
			  do i = 1, nb
				 
				 ib = i + 1
				 if(ib.gt.3) ib = 1
				 
				 pxab = pxi(i) - pxi(ib)
				 pyab = pyi(i) - pyi(ib)
				 pzab = pzi(i) - pzi(ib)
				 
				 pabsq = pxab**2 + pyab**2 + pzab**2
				 rabsq = dsqw(i)
				 diffsq = rabsq - pabsq
				 
				 if(abs(diffsq).gt.(rabsq * tol2)) then
					
					rxab = rxi(i) - rxi(ib)
					ryab = ryi(i) - ryi(ib)
					rzab = rzi(i) - rzi(ib)
					
					rpab = rxab * pxab + ryab * pyab + rzab * pzab
					
					rma = 1.0D0 / emw(i)
					rmb = 1.0D0 / emw(ib)
					gab = diffsq / (2.0D0 * (rma + rmb) * rpab)
					dx = rxab * gab
					dy = ryab * gab
					dz = rzab * gab
					
					pxi(i)  = pxi(i) + rma * dx
					pyi(i)  = pyi(i) + rma * dy
					pzi(i)  = pzi(i) + rma * dz
					pxi(ib) = pxi(ib) - rmb * dx
					pyi(ib) = pyi(ib) - rmb * dy
					pzi(ib) = pzi(ib) - rmb * dz
					
					vxi(i)  = vxi(i) + rma * dx
					vyi(i)  = vyi(i) + rma * dy
					vzi(i)  = vzi(i) + rma * dz
					vxi(ib) = vxi(ib) - rmb * dx
					vyi(ib) = vyi(ib) - rmb * dy
					vzi(ib) = vzi(ib) - rmb * dz
					
					ldone = .false.
					
				 endif
				 
			  enddo
			  
			  it = it + 1
			  goto 1000
			  
		   endif
		   
		   if(.not.ldone) then
			  
			  write(*,*) "iaq",iaq,"molecule ", j, "; Did not finish in move_w1_a"
			  open(10, file = "ERROR.txt")
			  write(10,*) "molecule ", j, "; Did not finish in move_w1_a"
			  close(10)
			  !stop
			  
		   endif
		   
		   do i = 1, 3
			  
			  iatom = (j - 1) * 3 + i
			  
			  rxw1(iatom, ib216) = pxi(i)
			  ryw1(iatom, ib216) = pyi(i)
			  rzw1(iatom, ib216) = pzi(i)
			  vxw1(iatom, ib4) = vxi(i)
			  vyw1(iatom, ib4) = vyi(i)
			  vzw1(iatom, ib4) = vzi(i)
			  
			enddo

		enddo
	
	enddo
	
	call nh_water1(w1v2)
	
	fxw1 = 0.0D0
	fyw1 = 0.0D0
	fzw1 = 0.0D0
	
	return
	end subroutine
