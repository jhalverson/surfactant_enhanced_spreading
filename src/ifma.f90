	integer*4 function igetIndex(iaq, ibq)

	implicit none
	integer*4 iaq, ibq

	igetIndex = iaq * (iaq + 1) + ibq + 1

	return
	end function


	complex*16 function Ylm(iaq, ibq, raq, rbq)

	implicit none
	integer*4 iaq, ibq, absm
	real*8 ri1, ri2, rifactorial
	real*8 raq, rbq, rlg, rLegendreP
      
	absm = abs(ibq)
	ri1 = rifactorial(iaq - absm)
	ri2 = rifactorial(iaq + absm)
	rlg = rLegendreP(iaq, absm, cos(raq))
	Ylm = sqrt(ri1 / ri2) * rlg * cmplx(cos(ibq * rbq), sin(ibq * rbq))
      
	return
	end function


	real*8 function rfctr(iaq, ibq)

	implicit none
	integer*4 iaq, ibq
	real*8 rifactorial

	rfctr = sqrt(rifactorial(iaq - abs(ibq)) / rifactorial(iaq + abs(ibq)))

	return
	end function


	real*8 function rAnm(iaq, ibq)

	implicit none
	integer iaq, ibq
	real*8 ri1, ri2, rifactorial

	ri1 = rifactorial(iaq - ibq)
	ri2 = rifactorial(iaq + ibq)
	rAnm = (-1)**iaq / sqrt(ri1 * ri2)
 
	return
	end function


	integer*4 function iIkm1(iaq, ibq)

	implicit none
	integer*4 iaq, ibq

	iIkm1 = (-1)**((abs(iaq) - abs(ibq) - abs(iaq - ibq)) / 2)

	return
	end function


	integer*4 function iIkm2(iaq, ibq)

	implicit none
	integer*4 iaq, ibq

	iIkm2 = (-1)**((abs(iaq - ibq) - abs(iaq) - abs(ibq)) / 2)

	return
	end function


	integer*4 function iIkm3(iaq, ibq)

	implicit none
	integer*4 iaq, ibq

	iIkm3 = (-1)**((abs(ibq) - abs(ibq - iaq) - abs(iaq)) / 2)

	return
	end function


	real*8 function rifactorial(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq

	rifactorial = rifact(iaq + 1)

	return
	end function
