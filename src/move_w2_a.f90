	subroutine move_w2_a(comm3Dq, iaq)
	
	implicit none
	include "mybox.cmns"
	include "mpif.h"
	
	integer*4 comm3Dq, iaq, ierr, ib512, ib64to512
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
	real*8 sum1, sum2
	real*8 sv1total, sv2total, w1v1total, w1v2total, w2v1total, w2v2total
	logical ldone
	
	tol = 0.0001D0
	maxit = 100
	
	nb = 3
	tol2 = 2.0D0 * tol
	
	escale = 1.0D0 - eta * dt / 2.0D0
	
	call nh_water2(w2v1)
	
	do ib4 = 1, 64
	
		ib512 = ib64to512(ib4)
		ictmolecule = ictw2(ib512)
		
		do j = 1, ictmolecule
		
		do i = 1, 3
			
			iatom = (j - 1) * 3 + i
			
			rxi(i) = rxw2(iatom, ib512)
			ryi(i) = ryw2(iatom, ib512)
			rzi(i) = rzw2(iatom, ib512)
			
			pxi(i) = rxw2(iatom, ib512) + vxw2(iatom, ib4) * escale + axw2(iatom, ib4)
			pyi(i) = ryw2(iatom, ib512) + vyw2(iatom, ib4) * escale + ayw2(iatom, ib4)
			pzi(i) = rzw2(iatom, ib512) + vzw2(iatom, ib4) * escale + azw2(iatom, ib4)
			
			vxi(i) = vxw2(iatom, ib4) * escale + axw2(iatom, ib4)
			vyi(i) = vyw2(iatom, ib4) * escale + ayw2(iatom, ib4)
			vzi(i) = vzw2(iatom, ib4) * escale + azw2(iatom, ib4)
			
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
			  
			  write(*,*) "iaq",iaq,"molecule ", j, "; Did not finish in move_w2_a"
			  open(10, file = "ERROR.txt")
			  write(10,*) "molecule ", j, "; Did not finish in move_w2_a"
			  close(10)
			  !stop
			  
		   endif
		   
		   do i = 1, 3
			  
			  iatom = (j - 1) * 3 + i
			  
			  rxw2(iatom, ib512) = pxi(i)
			  ryw2(iatom, ib512) = pyi(i)
			  rzw2(iatom, ib512) = pzi(i)
			  vxw2(iatom, ib4) = vxi(i)
			  vyw2(iatom, ib4) = vyi(i)
			  vzw2(iatom, ib4) = vzi(i)
			  
		   enddo
		   
		enddo
	
	enddo
	
	call nh_water2(w2v2)
	
	call mpi_allreduce(sv1, sv1total, 1, mpi_real8, mpi_sum, comm3Dq, ierr)
	call mpi_allreduce(sv2, sv2total, 1, mpi_real8, mpi_sum, comm3Dq, ierr)
	
	call mpi_allreduce(w1v1, w1v1total, 1, mpi_real8, mpi_sum, comm3Dq, ierr)
	call mpi_allreduce(w1v2, w1v2total, 1, mpi_real8, mpi_sum, comm3Dq, ierr)
	
	call mpi_allreduce(w2v1, w2v1total, 1, mpi_real8, mpi_sum, comm3Dq, ierr)
	call mpi_allreduce(w2v2, w2v2total, 1, mpi_real8, mpi_sum, comm3Dq, ierr)
	
	sum1 = (sv1total + w1v1total + w2v1total) / dt**2
	sum2 = (sv2total + w1v2total + w2v2total) / dt**2
	
	eta = eta + pnh1 * (sum1 + pnh2 + sum2)
	
	fxw2 = 0.0D0
	fyw2 = 0.0D0
	fzw2 = 0.0D0
	
	return
	end subroutine
