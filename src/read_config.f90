	subroutine read_config(iaq)
	
	implicit none
	include "mybox.cmns"
	
	integer*4 iaq, iend
	integer*4 ibox2, ibox64
	integer*4 ibO, ibG, ibH
	integer*4 ib64O, ib64G, ib64H
	integer*4 ib216, ib64to216, ib512, ib64to512, ib1000, ib64to1000
	integer*4 i, j, k
	
	real*8 xO, yO, zO, xG, yG, zG, xH, yH, zH
	real*8 rax(87532), rby(87532), rcz(87532)

	iend = 2 * (3 * iwater + 29 * isurfactant)
	
	! read atomic positions and velocities
	open(10, file = "data_cont2.dat")
	do i = 1, iend, 2

		! read and check position
		read(10,'(f8.3,f9.3,f9.3)') rax(i), rby(i), rcz(i)
		if(rax(i).lt.-rhsx.or.rax(i).gt.rhsx) write(*,*) "out of bounds"
		if(rby(i).lt.-rhsy.or.rby(i).gt.rhsy) write(*,*) "out of bounds"
		if(rcz(i).lt.-rhsz.or.rcz(i).gt.rhsz) write(*,*) "out of bounds"
		rax(i) = rax(i) + exp(-9.0D0)
		rby(i) = rby(i) + exp(-9.0D0)
		rcz(i) = rcz(i) + exp(-9.0D0)

		! read velocity
		read(10,'(f8.3,f9.3,f9.3)') rax(i + 1), rby(i + 1), rcz(i + 1)
		rax(i + 1) = (rax(i + 1) - 0.5D0) * 0.01D0
		rby(i + 1) = (rby(i + 1) - 0.5D0) * 0.01D0
		rcz(i + 1) = (rcz(i + 1) - 0.5D0) * 0.01D0

	enddo
	close(10)

	ictw1 = 0
	ictw2 = 0
	icts = 0

	! assign water molecules to the primary or split list
	do i = 1, 2 * (3 * iwater) - 5, 6

		xO = rax(i + 0)
		yO = rby(i + 0)
		zO = rcz(i + 0)

		xG = rax(i + 2)
		yG = rby(i + 2)
		zG = rcz(i + 2)

		xH = rax(i + 4)
		yH = rby(i + 4)
		zH = rcz(i + 4)

		ibO = ibox2(xO, yO, zO)
		ibG = ibox2(xG, yG, zG)
		ibH = ibox2(xH, yH, zH)

		if(ibO.eq.iaq + 1) then
		
			if(ibO.eq.ibG.and.ibO.eq.ibH) then
			
				ib64O = ibox64(xO, yO, zO)
				!if(ib64O.eq.64) write(*,*) xO, yO, zO
				ib64G = ibox64(xG, yG, zG)
				!if(ib64G.eq.64) write(*,*) xG, yG, zG
				ib64H = ibox64(xH, yH, zH)
				!if(ib64H.eq.64) write(*,*) xH, yH, zH

				if(ib64O.eq.ib64G.and.ib64O.eq.ib64H) then
				
					ib216 = ib64to216(ib64O)
					ictw1(ib216) = ictw1(ib216) + 1
					
					do j = 1, 6, 2
						
						k = (ictw1(ib216) - 1) * 3 + (j - 1) / 2 + 1

						rxw1(k, ib216) = rax(i + j - 1)
						ryw1(k, ib216) = rby(i + j - 1)
						rzw1(k, ib216) = rcz(i + j - 1)

						if(iaq.eq.1) then
						!write(*,*) rxw1(k, ib216), ryw1(k, ib216), rzw1(k, ib216)
						endif

						vxw1(k, ib64O) = rax(i + j)
						vyw1(k, ib64O) = rby(i + j)
						vzw1(k, ib64O) = rcz(i + j)

					enddo
					
				else
				
					ib512 = ib64to512(ib64O)
					ictw2(ib512) = ictw2(ib512) + 1

					do j = 1, 6, 2

						k = (ictw2(ib512) - 1) * 3 + (j - 1) / 2 + 1

						rxw2(k, ib512) = rax(i + j - 1)
						ryw2(k, ib512) = rby(i + j - 1)
						rzw2(k, ib512) = rcz(i + j - 1)

						vxw2(k, ib64O) = rax(i + j)
						vyw2(k, ib64O) = rby(i + j)
						vzw2(k, ib64O) = rcz(i + j)

					enddo
					
				endif
				
			else
			
				ib64O = ibox64(xO, yO, zO)
				ib512 = ib64to512(ib64O)
				ictw2(ib512) = ictw2(ib512) + 1

				do j = 1, 6, 2

					k = (ictw2(ib512) - 1) * 3 + (j - 1) / 2 + 1

					rxw2(k, ib512) = rax(i + j - 1)
					ryw2(k, ib512) = rby(i + j - 1)
					rzw2(k, ib512) = rcz(i + j - 1)

					vxw2(k, ib64O) = rax(i + j)
					vyw2(k, ib64O) = rby(i + j)
					vzw2(k, ib64O) = rcz(i + j)
					
				enddo
			
			endif
		
		endif
	
	enddo

	! assign surfactant molecules
	do i = 2 * (3 * iwater) + 27, iend - 31, 58

		xO = rax(i)
		yO = rby(i)
		zO = rcz(i)

		ibO = ibox2(xO, yO, zO)

		if(ibO.eq.iaq + 1) then

			ib64O = ibox64(xO, yO, zO)
			ib1000 = ib64to1000(ib64O)
			icts(ib1000) = icts(ib1000) + 1
			!if(i.eq.5607) write(*,*) iaq, ib64O, ib1000
			do j = 1, 58, 2

				k = (icts(ib1000) - 1) * 29 + (j - 1) / 2 + 1

				rxs(k, ib1000) = rax(i + j - 27)
				rys(k, ib1000) = rby(i + j - 27)
				rzs(k, ib1000) = rcz(i + j - 27)

				vxs(k, ib64O) = rax(i + j - 27 + 1)
				vys(k, ib64O) = rby(i + j - 27 + 1)
				vzs(k, ib64O) = rcz(i + j - 27 + 1)

			enddo

		endif

	enddo

	return
	end subroutine
