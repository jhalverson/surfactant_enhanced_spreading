	subroutine move_surfactant_b(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, ib1000, ib64to1000
	integer*4 nb, i, j, ictmolecule
	integer*4 iatom, ib4, it, ib, jb, maxit
	real*8 tol, tol2
	real*8 vxab, vyab, vzab
	real*8 rxab, ryab, rzab, rvab, rma, rmb, gab
	real*8 dx, dy, dz
	real*8 rxi(29), ryi(29), rzi(29)
	real*8 vxi(29), vyi(29), vzi(29)
	real*8 escale
	logical ldone
	
	tol = 0.0001D0
	maxit = 100
	
	nb = 28
	
	escale = 2.0D0 / (2.0D0 + eta * dt)
	
	do ib4 = 1, 64
	
		ib1000 = ib64to1000(ib4)
		
		ictmolecule = icts(ib1000)
		do j = 1, ictmolecule
		   
		do i = 1, 29
			  
			iatom = 29 * (j - 1) + i
			  
			rxi(i) = rxs(iatom, ib1000)
			ryi(i) = rys(iatom, ib1000)
			rzi(i) = rzs(iatom, ib1000)
			  
			axs(iatom, ib4) = fxs(iatom, ib1000) * dtsq2 / ems(i)
			ays(iatom, ib4) = fys(iatom, ib1000) * dtsq2 / ems(i)
			azs(iatom, ib4) = fzs(iatom, ib1000) * dtsq2 / ems(i)
			  
			vxi(i) = escale * (vxs(iatom, ib4) + axs(iatom, ib4))
			vyi(i) = escale * (vys(iatom, ib4) + ays(iatom, ib4))
			vzi(i) = escale * (vzs(iatom, ib4) + azs(iatom, ib4))
			  
		enddo
		   
		it = 0
		ldone = .false.
		   
1000	   if((.not.ldone).and.(it.le.maxit)) then
			  
			  ldone = .true.
			  
			  do i = 1, nb
				 
				 ib = ibond(i)
				 jb = jbond(i)
				 
				 vxab = vxi(ib) - vxi(jb)
				 vyab = vyi(ib) - vyi(jb)
				 vzab = vzi(ib) - vzi(jb)
				 
				 rxab = rxi(ib) - rxi(jb)
				 ryab = ryi(ib) - ryi(jb)
				 rzab = rzi(ib) - rzi(jb)
				 
				 rvab = rxab * vxab + ryab * vyab + rzab * vzab
				 rvab = rvab / dt
				 rma = 1.0D0 / ems(ib)
				 rmb = 1.0D0 / ems(jb)
				 gab = -rvab / ((rma + rmb) * dsqs(i))
				 
				 if(abs(gab).gt.tol) then
					
					dx = rxab * gab
					dy = ryab * gab
					dz = rzab * gab
					
					vxi(ib) = vxi(ib) + rma * dx * dt
					vyi(ib) = vyi(ib) + rma * dy * dt
					vzi(ib) = vzi(ib) + rma * dz * dt
					
					vxi(jb) = vxi(jb) - rmb * dx * dt
					vyi(jb) = vyi(jb) - rmb * dy * dt
					vzi(jb) = vzi(jb) - rmb * dz * dt
					
					ldone = .false.
					
				 endif
				 
			  enddo
			  
			  it = it + 1
			  goto 1000
			  
		   endif
		   
		   if(.not.ldone) then
			  
			  write(*,*) "iaq", iaq, "Too many iterations in move_surfactant_b: molecule", j
			  open(10, file = "ERROR.txt")
			  write(10,*) j, "Too many iterations in move_surfactant_b: molecule", j
			  close(10)
			  !stop
			  
		   endif
		   
		   do i = 1, 29
			  
			  iatom = 29 * (j - 1) + i
			  vxs(iatom, ib4) = vxi(i)
			  vys(iatom, ib4) = vyi(i)
			  vzs(iatom, ib4) = vzi(i)
			  
		   enddo
		   
		enddo
	
	enddo
	
	return
	end subroutine
