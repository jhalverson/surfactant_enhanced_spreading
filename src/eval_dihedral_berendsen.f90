	subroutine eval_dihedral_berendsen(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, ib4, i, j, ib1000, ib64to1000, ictmolecule
	integer*4 idv(30), jdv(30), kdv(30), ldv(30)
	integer*4 iv, jv, kv, lv
	real*8 rxji, ryji, rzji
	real*8 rxkj, rykj, rzkj
	real*8 rxnk, rynk, rznk
	real*8 cpjikjx, cpjikjy, cpjikjz
	real*8 cpkjnkx, cpkjnky, cpkjnkz
	real*8 cpjikjn, cpkjnkn, dp, arg
	real*8 thijkn, costh, sinth, r1
	real*8 fix, fiy, fiz, fix2, fiy2, fiz2, fxi, fyi, fzi
	real*8 fjx, fjy, fjz, fjx1, fjy1, fjz1, fjx2, fjy2, fjz2
	real*8 fjx3, fjy3, fjz3, fjx4, fjy4, fjz4, fjx5, fjy5, fjz5, fxj, fyj, fzj
	real*8 fkx, fky, fkz, fkx2, fky2, fkz2, fkx3, fky3, fkz3
	real*8 fkx4, fky4, fkz4, fkx5, fky5, fkz5, fkx6, fky6, fkz6, fxk, fyk, fzk
	real*8 fnx, fny, fnz, fnx2, fny2, fnz2, fxn, fyn, fzn
	
	do i = 1, 16
	   
	   idv(i) = i
	   jdv(i) = i + 1
	   kdv(i) = i + 2
	   ldv(i) = i + 3
	   
	enddo
	
	idv(17) = 16
	jdv(17) = 17
	kdv(17) = 18
	ldv(17) = 20
	
	idv(18) = 16
	jdv(18) = 17
	kdv(18) = 18
	ldv(18) = 25
	
	idv(19) = 17
	jdv(19) = 18
	kdv(19) = 20
	ldv(19) = 21
	
	idv(20) = 17
	jdv(20) = 18
	kdv(20) = 25
	ldv(20) = 26
	
	idv(21) = 19
	jdv(21) = 18
	kdv(21) = 20
	ldv(21) = 21
	
	idv(22) = 19
	jdv(22) = 18
	kdv(22) = 25
	ldv(22) = 26
	
	idv(23) = 18
	jdv(23) = 20
	kdv(23) = 21
	ldv(23) = 22
	
	idv(24) = 18
	jdv(24) = 20
	kdv(24) = 21
	ldv(24) = 23
	
	idv(25) = 18
	jdv(25) = 20
	kdv(25) = 21
	ldv(25) = 24
	
	idv(26) = 18
	jdv(26) = 25
	kdv(26) = 26
	ldv(26) = 27
	
	idv(27) = 18
	jdv(27) = 25
	kdv(27) = 26
	ldv(27) = 28
	
	idv(28) = 18
	jdv(28) = 25
	kdv(28) = 26
	ldv(28) = 29
	
	idv(29) = 20
	jdv(29) = 18
	kdv(29) = 25
	ldv(29) = 26
	
	idv(30) = 21
	jdv(30) = 20
	kdv(30) = 18
	ldv(30) = 25
	
	do ib4 = 1, 64
	
		ib1000 = ib64to1000(ib4)
		ictmolecule = icts(ib1000)
		
		do j = 1, ictmolecule
		
			do i = 19, 30
			
				iv = (j - 1) * 29 + idv(i)
				jv = (j - 1) * 29 + jdv(i)
				kv = (j - 1) * 29 + kdv(i)
				lv = (j - 1) * 29 + ldv(i)
				
				rxji = rxs(jv, ib1000) - rxs(iv, ib1000)
				ryji = rys(jv, ib1000) - rys(iv, ib1000)
				rzji = rzs(jv, ib1000) - rzs(iv, ib1000)
				
				rxkj = rxs(kv, ib1000) - rxs(jv, ib1000)
				rykj = rys(kv, ib1000) - rys(jv, ib1000)
				rzkj = rzs(kv, ib1000) - rzs(jv, ib1000)
				
				rxnk = rxs(lv, ib1000) - rxs(kv, ib1000)
				rynk = rys(lv, ib1000) - rys(kv, ib1000)
				rznk = rzs(lv, ib1000) - rzs(kv, ib1000)
				
				cpjikjx = rxji * rykj - ryji * rxkj
				cpjikjy = rxji * rzkj - rzji * rxkj
				cpjikjz = ryji * rzkj - rzji * rykj
	      
				cpkjnkx = rxkj * rynk - rykj * rxnk
				cpkjnky = rxkj * rznk - rzkj * rxnk
				cpkjnkz = rykj * rznk - rzkj * rynk
	      
				cpjikjn = dsqrt(cpjikjx**2 + cpjikjy**2 + cpjikjz**2)
				cpkjnkn = dsqrt(cpkjnkx**2 + cpkjnky**2 + cpkjnkz**2)
	      
				dp = cpjikjx * cpkjnkx + cpjikjy * cpkjnky + cpjikjz * cpkjnkz
				arg = dp / cpjikjn / cpkjnkn
	      
				if(arg.le.-1.0D0) arg = -0.999999D0
				if(arg.ge.1.0D0) arg = 0.999999D0
	      
				thijkn = acos(arg)
				costh = cos(thijkn)
				sinth = sin(thijkn)
	      
				r1 = -3.0D0 * dva(i) * sin(3.0D0 * thijkn) / sinth
	      
				! particle i
	      
				fix = rxkj * (-rykj * rynk - rzkj * rznk) + rxnk * (rykj * rykj + rzkj * rzkj)
				fiy = rykj * (-rxkj * rxnk - rzkj * rznk) + rynk * (rxkj * rxkj + rzkj * rzkj)
				fiz = rzkj * (-rxkj * rxnk - rykj * rynk) + rznk * (rxkj * rxkj + rykj * rykj)
	      
				fix = fix / cpjikjn / cpkjnkn
				fiy = fiy / cpjikjn / cpkjnkn
				fiz = fiz / cpjikjn / cpkjnkn
	      
				fix2 = 2.0D0 * rxji * (-rykj * rykj - rzkj * rzkj) + 2.0D0 * rxkj * (ryji * rykj + rzji * rzkj)
				fiy2 = 2.0D0 * ryji * (-rxkj * rxkj - rzkj * rzkj) + 2.0D0 * rykj * (rxji * rxkj + rzji * rzkj)
				fiz2 = 2.0D0 * rzji * (-rxkj * rxkj - rykj * rykj) + 2.0D0 * rzkj * (rxji * rxkj + ryji * rykj)
	      
				fix2 = fix2 / cpjikjn**2
				fiy2 = fiy2 / cpjikjn**2
				fiz2 = fiz2 / cpjikjn**2
	      
				fxi = r1 * (fix - 0.5D0 * costh * fix2)
				fyi = r1 * (fiy - 0.5D0 * costh * fiy2)
				fzi = r1 * (fiz - 0.5D0 * costh * fiz2)
	      
				! particle j
	      
				fjx = rxji * (-rykj * rynk - rzkj * rznk) + rxkj * (rykj * rynk + rzkj * rznk)
				fjy = ryji * (-rxkj * rxnk - rzkj * rznk) + rykj * (rxkj * rxnk + rzkj * rznk)
				fjz = rzji * (-rxkj * rxnk - rykj * rynk) + rzkj * (rxkj * rxnk + rykj * rynk)
				fjx1 = rxnk * (-ryji * rykj - rzji * rzkj - rykj * rykj - rzkj * rzkj)
				fjy1 = rynk * (-rxji * rxkj - rzji * rzkj - rxkj * rxkj - rzkj * rzkj)
				fjz1 = rznk * (-rxji * rxkj - ryji * rykj - rxkj * rxkj - rykj * rykj)
				fjx2 = 2.0D0 * rxkj * (ryji * rynk + rzji * rznk)
				fjy2 = 2.0D0 * rykj * (rxji * rxnk + rzji * rznk)
				fjz2 = 2.0D0 * rzkj * (rxji * rxnk + ryji * rynk)
	      
				fjx = (fjx + fjx1 + fjx2) / cpjikjn / cpkjnkn
				fjy = (fjy + fjy1 + fjy2) / cpjikjn / cpkjnkn
				fjz = (fjz + fjz1 + fjz2) / cpjikjn / cpkjnkn
	      
				fjx3 = 2.0D0 * rxji * (rykj * rykj + rzkj * rzkj + ryji * rykj + rzji * rzkj)
				fjy3 = 2.0D0 * ryji * (rxkj * rxkj + rzkj * rzkj + rxji * rxkj + rzji * rzkj)
				fjz3 = 2.0D0 * rzji * (rxkj * rxkj + rykj * rykj + rxji * rxkj + ryji * rykj)
				fjx4 = 2.0D0 * rxkj * (-ryji * ryji - rzji * rzji - ryji * rykj - rzji * rzkj)
				fjy4 = 2.0D0 * rykj * (-rxji * rxji - rzji * rzji - rxji * rxkj - rzji * rzkj)
				fjz4 = 2.0D0 * rzkj * (-rxji * rxji - ryji * ryji - rxji * rxkj - ryji * rykj)
	      
				fjx3 = (fjx3 + fjx4) / cpjikjn**2
				fjy3 = (fjy3 + fjy4) / cpjikjn**2
				fjz3 = (fjz3 + fjz4) / cpjikjn**2
	      
				fjx5 = 2.0D0 * rxnk * (rykj * rynk + rzkj * rznk) + 2.0D0 * rxkj * (-rynk * rynk - rznk * rznk)
				fjy5 = 2.0D0 * rynk * (rxkj * rxnk + rzkj * rznk) + 2.0D0 * rykj * (-rxnk * rxnk - rznk * rznk)
				fjz5 = 2.0D0 * rznk * (rxkj * rxnk + rykj * rynk) + 2.0D0 * rzkj * (-rxnk * rxnk - rynk * rynk)
	      
				fjx5 = fjx5 / cpkjnkn**2
				fjy5 = fjy5 / cpkjnkn**2
				fjz5 = fjz5 / cpkjnkn**2
	      
				fxj = r1 * (fjx - 0.5D0 * costh * (fjx3 + fjx5))
				fyj = r1 * (fjy - 0.5D0 * costh * (fjy3 + fjy5))
				fzj = r1 * (fjz - 0.5D0 * costh * (fjz3 + fjz5))
	      
				! particle k
	      
				fkx = rxji * (rykj * rykj + rzkj * rzkj + rykj * rynk + rzkj * rznk)
				fky = ryji * (rxkj * rxkj + rzkj * rzkj + rxkj * rxnk + rzkj * rznk)
				fkz = rzji * (rxkj * rxkj + rykj * rykj + rxkj * rxnk + rykj * rynk)
				fkx2 = rxkj * (-ryji * rykj - rzji * rzkj)
				fky2 = rykj * (-rxji * rxkj - rzji * rzkj)
				fkz2 = rzkj * (-rxji * rxkj - ryji * rykj)
				fkx3 = rxnk * (ryji * rykj + rzji * rzkj) + 2.0D0 * rxkj * (-ryji * rynk - rzji * rznk)
				fky3 = rynk * (rxji * rxkj + rzji * rzkj) + 2.0D0 * rykj * (-rxji * rxnk - rzji * rznk)
				fkz3 = rznk * (rxji * rxkj + ryji * rykj) + 2.0D0 * rzkj * (-rxji * rxnk - ryji * rynk)
	      
				fkx = (fkx + fkx2 + fkx3) / cpjikjn / cpkjnkn
				fky = (fky + fky2 + fky3) / cpjikjn / cpkjnkn
				fkz = (fkz + fkz2 + fkz3) / cpjikjn / cpkjnkn
	      
				fkx4 = 2.0D0 * rxji * (-ryji * rykj - rzji * rzkj) + 2.0D0 * rxkj * (ryji * ryji + rzji * rzji)
				fky4 = 2.0D0 * ryji * (-rxji * rxkj - rzji * rzkj) + 2.0D0 * rykj * (rxji * rxji + rzji * rzji)
				fkz4 = 2.0D0 * rzji * (-rxji * rxkj - ryji * rykj) + 2.0D0 * rzkj * (rxji * rxji + ryji * ryji)
	      
				fkx4 = fkx4 / cpjikjn**2
				fky4 = fky4 / cpjikjn**2
				fkz4 = fkz4 / cpjikjn**2
	      
				fkx5 = 2.0D0 * rxnk * (-rykj * rykj - rzkj * rzkj - rykj * rynk - rzkj * rznk)
				fky5 = 2.0D0 * rynk * (-rxkj * rxkj - rzkj * rzkj - rxkj * rxnk - rzkj * rznk)
				fkz5 = 2.0D0 * rznk * (-rxkj * rxkj - rykj * rykj - rxkj * rxnk - rykj * rynk)
				fkx6 = 2.0D0 * rxkj * (rynk * rynk + rznk * rznk + rykj * rynk + rzkj * rznk)
				fky6 = 2.0D0 * rykj * (rxnk * rxnk + rznk * rznk + rxkj * rxnk + rzkj * rznk)
				fkz6 = 2.0D0 * rzkj * (rxnk * rxnk + rynk * rynk + rxkj * rxnk + rykj * rynk)
	      
				fkx5 = (fkx5 + fkx6) / cpkjnkn**2
				fky5 = (fky5 + fky6) / cpkjnkn**2
				fkz5 = (fkz5 + fkz6) / cpkjnkn**2
	      
				fxk = r1 * (fkx - 0.5D0 * costh * (fkx4 + fkx5))
				fyk = r1 * (fky - 0.5D0 * costh * (fky4 + fky5))
				fzk = r1 * (fkz - 0.5D0 * costh * (fkz4 + fkz5))
	      
				! particle n
	      
				fnx = rxji * (-rykj * rykj - rzkj * rzkj) + rxkj * (ryji * rykj + rzji * rzkj)
				fny = ryji * (-rxkj * rxkj - rzkj * rzkj) + rykj * (rxji * rxkj + rzji * rzkj)
				fnz = rzji * (-rxkj * rxkj - rykj * rykj) + rzkj * (rxji * rxkj + ryji * rykj)
	      
				fnx = fnx / cpjikjn / cpkjnkn
				fny = fny / cpjikjn / cpkjnkn
				fnz = fnz / cpjikjn / cpkjnkn
	      
				fnx2 = 2.0D0 * rxnk * (rykj * rykj + rzkj * rzkj) + 2.0D0 * rxkj * (-rykj * rynk - rzkj * rznk)
				fny2 = 2.0D0 * rynk * (rxkj * rxkj + rzkj * rzkj) + 2.0D0 * rykj * (-rxkj * rxnk - rzkj * rznk)
				fnz2 = 2.0D0 * rznk * (rxkj * rxkj + rykj * rykj) + 2.0D0 * rzkj * (-rxkj * rxnk - rykj * rynk)
	      
				fnx2 = fnx2 / cpkjnkn**2
				fny2 = fny2 / cpkjnkn**2
				fnz2 = fnz2 / cpkjnkn**2
	      
				fxn = r1 * (fnx - 0.5D0 * costh * fnx2)
				fyn = r1 * (fny - 0.5D0 * costh * fny2)
				fzn = r1 * (fnz - 0.5D0 * costh * fnz2)
	      
				! forces i, j, k, n
				
				fxs(iv, ib1000) = fxs(iv, ib1000) + fxi
				fys(iv, ib1000) = fys(iv, ib1000) + fyi
				fzs(iv, ib1000) = fzs(iv, ib1000) + fzi
				
				fxs(jv, ib1000) = fxs(jv, ib1000) + fxj
				fys(jv, ib1000) = fys(jv, ib1000) + fyj
				fzs(jv, ib1000) = fzs(jv, ib1000) + fzj
				
				fxs(kv, ib1000) = fxs(kv, ib1000) + fxk
				fys(kv, ib1000) = fys(kv, ib1000) + fyk
				fzs(kv, ib1000) = fzs(kv, ib1000) + fzk
				
				fxs(lv, ib1000) = fxs(lv, ib1000) + fxn
				fys(lv, ib1000) = fys(lv, ib1000) + fyn
				fzs(lv, ib1000) = fzs(lv, ib1000) + fzn
				
			enddo
			
		enddo
		
	enddo
	
	return
	end subroutine
