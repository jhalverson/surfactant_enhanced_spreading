	subroutine eval_valence(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, i, j, ib4, ib1000, ib64to1000, ictmolecule
	integer*4 iang(36), jang(36), kang(36)
	integer*4 ig, jg, kg
	real*8 rxji, ryji, rzji, rxki, ryki, rzki
	real*8 dp, rjin, rkin, thjik, costh, r1
	real*8 fjx, fjy, fjz, fkx, fky, fkz
	
	do i = 1, 17
	   
	   iang(i) = i
	   jang(i) = i + 1
	   kang(i) = i + 2
	   
	enddo
	
	iang(18) = 17
	jang(18) = 18
	kang(18) = 20
	
	iang(19) = 17
	jang(19) = 18
	kang(19) = 25
	
	iang(20) = 19
	jang(20) = 18
	kang(20) = 20
	
	iang(21) = 19
	jang(21) = 18
	kang(21) = 25
	
	iang(22) = 18
	jang(22) = 25
	kang(22) = 26
	
	iang(23) = 18
	jang(23) = 20
	kang(23) = 21
	
	iang(24) = 20
	jang(24) = 21
	kang(24) = 22
	
	iang(25) = 20
	jang(25) = 21
	kang(25) = 23
	
	iang(26) = 20
	jang(26) = 21
	kang(26) = 24
	
	iang(27) = 22
	jang(27) = 21
	kang(27) = 24
	
	iang(28) = 22
	jang(28) = 21
	kang(28) = 23
	
	iang(29) = 23
	jang(29) = 21
	kang(29) = 24
	
	iang(30) = 25
	jang(30) = 26
	kang(30) = 27
	
	iang(31) = 25
	jang(31) = 26
	kang(31) = 28
	
	iang(32) = 25
	jang(32) = 26
	kang(32) = 29
	
	iang(33) = 27
	jang(33) = 26
	kang(33) = 28
	
	iang(34) = 27
	jang(34) = 26
	kang(34) = 29
	
	iang(35) = 28
	jang(35) = 26
	kang(35) = 29
	
	iang(36) = 20
	jang(36) = 18
	kang(36) = 25
	
	do ib4 = 1, 64
	
		ib1000 = ib64to1000(ib4)
		ictmolecule = icts(ib1000)
		
		do j = 1, ictmolecule

			do i = 1, 36

				ig = (j - 1) * 29 + iang(i)
				jg = (j - 1) * 29 + jang(i)
				kg = (j - 1) * 29 + kang(i)

				rxji = rxs(ig, ib1000) - rxs(jg, ib1000)
				ryji = rys(ig, ib1000) - rys(jg, ib1000)
				rzji = rzs(ig, ib1000) - rzs(jg, ib1000)

				rxki = rxs(kg, ib1000) - rxs(jg, ib1000)
				ryki = rys(kg, ib1000) - rys(jg, ib1000)
				rzki = rzs(kg, ib1000) - rzs(jg, ib1000)
				
				dp = rxji * rxki + ryji * ryki + rzji * rzki
				rjin = dsqrt(rxji**2 + ryji**2 + rzji**2)
				rkin = dsqrt(rxki**2 + ryki**2 + rzki**2)
				  
				thjik = acos(dp / rjin / rkin)
				costh = cos(thjik)
				  
				r1 = rkth(i) * (thjik - th(i)) / sin(thjik)
				  
				fjx = rxki / rjin / rkin - costh * rxji / rjin**2
				fjy = ryki / rjin / rkin - costh * ryji / rjin**2
				fjz = rzki / rjin / rkin - costh * rzji / rjin**2
				  
				fkx = rxji / rjin / rkin - costh * rxki / rkin**2
				fky = ryji / rjin / rkin - costh * ryki / rkin**2
				fkz = rzji / rjin / rkin - costh * rzki / rkin**2
				
				fxs(ig, ib1000) = fxs(ig, ib1000) + r1 * fjx
				fys(ig, ib1000) = fys(ig, ib1000) + r1 * fjy
				fzs(ig, ib1000) = fzs(ig, ib1000) + r1 * fjz
				
				fxs(jg, ib1000) = fxs(jg, ib1000) - r1 * (fjx + fkx)
				fys(jg, ib1000) = fys(jg, ib1000) - r1 * (fjy + fky)
				fzs(jg, ib1000) = fzs(jg, ib1000) - r1 * (fjz + fkz)
				
				fxs(kg, ib1000) = fxs(kg, ib1000) + r1 * fkx
				fys(kg, ib1000) = fys(kg, ib1000) + r1 * fky
				fzs(kg, ib1000) = fzs(kg, ib1000) + r1 * fkz
				
				if(iaq.eq.14.and.ib1000.eq.556.and.j.eq.1) then
	!write(*,'(f10.3,f10.3,f10.3,f10.3,f10.3,f10.3)') fxs(ig, ib1000),fys(ig, ib1000),fzs(ig, ib1000),rxs(ig, ib1000),rys(ig,ib1000),rzs(ig,ib1000)
	!write(*,'(f10.3,f10.3,f10.3,f10.3,f10.3,f10.3)') fxs(jg, ib1000),fys(jg, ib1000),fzs(jg, ib1000),rxs(jg, ib1000),rys(jg,ib1000),rzs(jg,ib1000)
	!write(*,'(f10.3,f10.3,f10.3,f10.3,f10.3,f10.3)') fxs(kg, ib1000),fys(kg, ib1000),fzs(kg, ib1000),rxs(kg, ib1000),rys(kg,ib1000),rzs(kg,ib1000)
	!write(*,*) rkth(i)
				endif
				
			enddo
			
		enddo
		
	enddo
	
	return
	end subroutine
