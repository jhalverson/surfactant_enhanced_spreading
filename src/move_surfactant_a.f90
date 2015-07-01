	subroutine move_surfactant_a(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, ib1000, ib64to1000
	integer*4 nb, i, j, ictmolecule
	integer*4 iatom, ib4, it, ib, jb, maxit
	real*8 tol, tol2
	real*8 pxab, pyab, pzab, pabsq, rabsq, diffsq
	real*8 rxab, ryab, rzab, rpab, rma, rmb, gab
	real*8 dx, dy, dz
	real*8 rxi(29), ryi(29), rzi(29)
	real*8 pxi(29), pyi(29), pzi(29)
	real*8 vxi(29), vyi(29), vzi(29)
	real*8 escale
	logical ldone
	
	tol = 0.0001D0
	maxit = 100
	
	nb = 28
	tol2 = 2.0D0 * tol
	
	escale = 1.0D0 - eta * dt / 2.0D0
	
	call nh_surfactant(sv1)
	
	do ib4 = 1, 64
	
		ib1000 = ib64to1000(ib4)
		ictmolecule = icts(ib1000)
		
		do j = 1, ictmolecule
		
			do i = 1, 29
			
			iatom = (j - 1) * 29 + i
			
			rxi(i) = rxs(iatom, ib1000)
			ryi(i) = rys(iatom, ib1000)
			rzi(i) = rzs(iatom, ib1000)
			
			pxi(i) = rxs(iatom, ib1000) + vxs(iatom, ib4) * escale + axs(iatom, ib4)
			pyi(i) = rys(iatom, ib1000) + vys(iatom, ib4) * escale + ays(iatom, ib4)
			pzi(i) = rzs(iatom, ib1000) + vzs(iatom, ib4) * escale + azs(iatom, ib4)
			
			vxi(i) = vxs(iatom, ib4) * escale + axs(iatom, ib4)
			vyi(i) = vys(iatom, ib4) * escale + ays(iatom, ib4)
			vzi(i) = vzs(iatom, ib4) * escale + azs(iatom, ib4)
			
		enddo
		   
		it = 0
		ldone = .false.
		   
1000	   if((.not.ldone).and.(it.le.maxit)) then
		   ldone = .true.

			  do i = 1, nb
				 
				 ib = ibond(i)
				 jb = jbond(i)
				 
				 pxab = pxi(ib) - pxi(jb)
				 pyab = pyi(ib) - pyi(jb)
				 pzab = pzi(ib) - pzi(jb)
				 
				 pabsq = pxab**2 + pyab**2 + pzab**2
				 rabsq = dsqs(i)
				 diffsq = rabsq - pabsq
				 
				 if(abs(diffsq).gt.(rabsq * tol2)) then
					
					rxab = rxi(ib) - rxi(jb)
					ryab = ryi(ib) - ryi(jb)
					rzab = rzi(ib) - rzi(jb)
					
					rpab = rxab * pxab + ryab * pyab + rzab * pzab
					
					rma = 1.0D0 / ems(ib)
					rmb = 1.0D0 / ems(jb)
					gab = diffsq / (2.0D0 * (rma + rmb) * rpab)
					dx = rxab * gab
					dy = ryab * gab
					dz = rzab * gab
					
					pxi(ib) = pxi(ib) + rma * dx
					pyi(ib) = pyi(ib) + rma * dy
					pzi(ib) = pzi(ib) + rma * dz
					pxi(jb) = pxi(jb) - rmb * dx
					pyi(jb) = pyi(jb) - rmb * dy
					pzi(jb) = pzi(jb) - rmb * dz
					
					vxi(ib) = vxi(ib) + rma * dx
					vyi(ib) = vyi(ib) + rma * dy
					vzi(ib) = vzi(ib) + rma * dz
					vxi(jb) = vxi(jb) - rmb * dx
					vyi(jb) = vyi(jb) - rmb * dy
					vzi(jb) = vzi(jb) - rmb * dz

					ldone = .false.

				 endif

			  enddo

			  it = it + 1
			  goto 1000

		   endif

		   if(.not.ldone) then

			  write(*,*) "iaq",iaq, "molecule ", j, "; Did not finish in move_surfactant_a"
			  open(10, file = "ERROR.txt")
			  write(10,*) "molecule ", j, "; Did not finish in move_surfactant_a"
			  close(10)
			  !stop

		   endif

		   do i = 1, 29

			  iatom = 29 * (j - 1) + i

			  rxs(iatom, ib1000) = pxi(i)
			  rys(iatom, ib1000) = pyi(i)
			  rzs(iatom, ib1000) = pzi(i)
			
			  vxs(iatom, ib4) = vxi(i)
			  vys(iatom, ib4) = vyi(i)
			  vzs(iatom, ib4) = vzi(i)

		   enddo

		enddo

	enddo
	
	call nh_surfactant(sv2)
	
	fxs = 0.0D0
	fys = 0.0D0
	fzs = 0.0D0
	
	return
	end subroutine
