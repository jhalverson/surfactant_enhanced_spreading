	program make_graphite
	
	implicit none
	integer*4 np, i, j, ix, iy, ins, iatom
	parameter (np = 1000000)
	real*8 x(np), y(np), z(np)
	real*8 pi, sigma, ccbond, xshift, yshift, zmin
	real*8 xcarbon, ycarbon, zcarbon, yup
	real*8 xmin, xmax, ymin, ymax
	
	character(20) flnm
	character(10) substr
	character(32) numform
	character(10) fldst
	
	numform = '(a6,i5,a5,a4,i6,f12.3,f8.3,f8.3)'

	pi = 4.0D0 * atan(1.0D0)
	sigma = 3.166D0
	ccbond = 1.421D0 / sigma
	xshift = cos(60.0D0 * pi / 180.0D0) * ccbond
	yshift = sin(60.0D0 * pi / 180.0D0) * ccbond

	! ix should be an odd integer
	ix = 81
	iy = sin(60.0D0 * pi / 180.0D0) * ix
	ins = 2 * ix * iy + 2 * (iy - 1)

	xcarbon = 0.0D0
	ycarbon = 0.0D0
	zcarbon = -16.0D0

	iatom = 0
	do i = 1, ix

		do j = 1, iy

			yup = 0.0D0
			if(mod(i, 2).eq.0) yup = yshift

			iatom = iatom + 1
			x(iatom) = i * (ccbond + xshift)
			y(iatom) = j * 2.0D0 * yshift + yup
			z(iatom) = zcarbon

			iatom = iatom + 1
			x(iatom) = x(iatom - 1) + ccbond
			y(iatom) = y(iatom - 1)
			z(iatom) = zcarbon
			
		enddo

	enddo

	do j = 1, iy - 1

		iatom = iatom + 1
		x(iatom) = ccbond
		y(iatom) = j * 2.0D0 * yshift + yshift
		z(iatom) = zcarbon

	enddo		

	do j = 1, iy - 1

		iatom = iatom + 1
		x(iatom) = (ix + 1) * (ccbond + xshift)
		y(iatom) = (j + 1) * 2.0D0 * yshift - yshift
		z(iatom) = zcarbon

	enddo

	do j = 1, ins

		iatom = iatom + 1
		x(iatom) = x(iatom - ins) + ccbond
		y(iatom) = y(iatom - ins)
		z(iatom) = z(iatom - ins) - 3.41D0 / sigma

	enddo
	
	xmin = x(1)
	xmax = x(1)
	ymin = y(1)
	ymax = y(1)
	
	do j = 2, iatom
	
		if(x(j).gt.xmax) xmax = x(j)
		if(x(j).lt.xmin) xmin = x(j)
		if(y(j).gt.ymax) ymax = y(j)
		if(y(j).lt.ymin) ymin = y(j)
	
	enddo
	
	do j = 1, iatom
	
		x(j) = x(j) - (xmax + xmin) / 2.0D0 + xcarbon
		y(j) = y(j) - (ymax + ymin) / 2.0D0 + ycarbon
		
	enddo
	
	xmin = x(1)
	xmax = x(1)
	ymin = y(1)
	ymax = y(1)
	
	do j = 2, iatom
	
		if(x(j).gt.xmax) xmax = x(j)
		if(x(j).lt.xmin) xmin = x(j)
		if(y(j).gt.ymax) ymax = y(j)
		if(y(j).lt.ymin) ymin = y(j)
	
	enddo
	
	write(*,*) "atoms = ", iatom
	write(*,*) "x-side =", (xmax - xmin) * sigma, "Angstrom"
	write(*,*) "y-side =", (ymax - ymin) * sigma, "Angstrom"
	write(*,*) "z-side =", 3.41D0, "Angstrom"
	write(*,*) "x-center =", 0.5D0 * (xmax + xmin) * sigma, "Angstrom"
	write(*,*) "y-center =", 0.5D0 * (ymax + ymin) * sigma, "Angstrom"
	write(*,*) "z-center =", zcarbon, "Angstrom"
	
	open(10, file = "graphite.dat")
	do j = 1, iatom
	write(10,'(f10.4,f11.4,f11.4)') x(j), y(j), z(j)
	enddo
	close(10)
	
	open(10, file = "graphite.pdb")
	do j = 1, iatom
	write(10, numform) "HETATM", j, "C","GPH", j, x(j)*sigma, y(j)*sigma, z(j)*sigma
	enddo
	close(10)

	end program
