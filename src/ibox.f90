	integer*4 function ib64to1000(iaq)

	implicit none
	integer*4 iaq, k, m, iaqr

	k = (iaq - 1) / 16
	m = (iaq - 16 * k - 1) / 4
	iaqr = iaq - 16 * k - 4 * m

	ib64to1000 = 333 + 100 * k + 10 * m + iaqr

	return
	end function


	integer*4 function ib64to512(iaq)

	implicit none
	integer*4 iaq, k, m, iaqr

	k = (iaq - 1) / 16
	m = (iaq - 16 * k - 1) / 4
	iaqr = iaq - 16 * k - 4 * m

	ib64to512 = 146 + 64 * k + 8 * m + iaqr

	return
	end function


	integer*4 function ib64to216(iaq)

	implicit none
	integer*4 iaq, k, m, iaqr

	k = (iaq - 1) / 16
	m = (iaq - 16 * k - 1) / 4
	iaqr = iaq - 16 * k - 4 * m

	ib64to216 = 43 + 36 * k + 6 * m + iaqr

	return
	end function


	integer*4 function ib8to216(iaq)

	implicit none
	integer*4 iaq, k, m, iaqr

	k = (iaq - 1) / 4
	m = (iaq - 4 * k - 1) / 2
	iaqr = iaq - 4 * k - 2 * m

	ib8to216 = 86 + 36 * k + 6 * m + iaqr

	return
	end function


	integer*4 function ib512to64(iaq)

	implicit none
	integer*4 iaq, k, m

	k = (iaq - 1) / 64
	m = (iaq - 64 * k - 1) / 8

	ib512to64 = iaq - 64 * k - 8 * m - 2 + 16 * (k - 2) + 4 * (m - 2)

	return
	end function


	logical function isnearest(raq, rbq, rcq, rdq, req, rfq)

	implicit none
	include "mybox.cmns"

	real*8 raq, rbq, rcq, rdq, req, rfq
	real*8 rs15
	logical boolx, booly, boolz

	rs15 = 1.5D0 * rs

	boolx = .false.
	booly = .false.
	boolz = .false.

	if(raq.lt.rdq + rs15.and.raq.ge.rdq - rs15) boolx = .true.
	if(rbq.lt.req + rs15.and.rbq.ge.req - rs15) booly = .true.
	if(rcq.lt.rfq + rs15.and.rcq.ge.rfq - rs15) boolz = .true.

	isnearest = .false.
	if(boolx.and.booly.and.boolz) isnearest = .true.

	return
	end function


	logical function isCentral(raq, rbq, rcq, rdq, req, rfq)

	implicit none
	include "mybox.cmns"

	real*8 raq, rbq, rcq, rdq, req, rfq
	real*8 rs2
	logical boolx, booly, boolz

	rs2 = rs / 2.0D0

	boolx = .false.
	booly = .false.
	boolz = .false.

	if(raq.lt.rdq + rs2.and.raq.ge.rdq - rs2) boolx = .true.
	if(rbq.lt.req + rs2.and.rbq.ge.req - rs2) booly = .true.
	if(rcq.lt.rfq + rs2.and.rcq.ge.rfq - rs2) boolz = .true.

	isCentral = .false.
	if(boolx.and.booly.and.boolz) isCentral = .true.

	return
	end function


	integer*4 function ibox2(raq, rbq, rcq)
	
	implicit none
	include "mybox.cmns"
	
	real*8 raq, rbq, rcq

	ibox2 = 1 + int((raq + rhsx) * ibdx2 / rsdx)
	ibox2 = ibox2 + int((rbq + rhsy) * ibdy2 / rsdy) * ibdx2
	ibox2 = ibox2 + int((rcq + rhsz) * ibdz2 / rsdz) * ibdx2 * ibdy2
	
	return
	end function


	integer*4 function ibox3(raq, rbq, rcq)
	
	implicit none
	include "mybox.cmns"
	
	real*8 raq, rbq, rcq
	
	ibox3 = 1 + int((raq + rhsx) * ibdx3 / rsdx)
	ibox3 = ibox3 + int((rbq + rhsy) * ibdy3 / rsdy) * ibdx3
	ibox3 = ibox3 + int((rcq + rhsz) * ibdz3 / rsdz) * ibdx3 * ibdy3

	return
	end function


	integer*4 function ibox4(raq, rbq, rcq)
	
	implicit none
	include "mybox.cmns"
	
	real*8 raq, rbq, rcq
	
	ibox4 = 1 + int((raq + rhsx) * ibdx4 / rsdx)
	ibox4 = ibox4 + int((rbq + rhsy) * ibdy4 / rsdy) * ibdx4
	ibox4 = ibox4 + int((rcq + rhsz) * ibdz4 / rsdz) * ibdx4 * ibdy4

	return
	end function


	integer*4 function ibox64(raq, rbq, rcq)

	implicit none
	include "mybox.cmns"
	real*8 raq, rbq, rcq
	real*8 rap, rbp, rcp

	rap = raq - rcen2xtmp + 2.0D0 * rs
	rbp = rbq - rcen2ytmp + 2.0D0 * rs
	rcp = rcq - rcen2ztmp + 2.0D0 * rs

	ibox64 = 1 + int(rap / rs) + int(rbp / rs) * 4 + int(rcp / rs) * 4 * 4

	return
	end function


	integer*4 function ibox64mod(raq, rbq, rcq, rdq, req, rfq)

	implicit none
	include "mybox.cmns"
	real*8 raq, rbq, rcq, rdq, req, rfq
	real*8 rap, rbp, rcp

	rap = raq - rdq + 2.0D0 * rs
	rbp = rbq - req + 2.0D0 * rs
	rcp = rcq - rfq + 2.0D0 * rs

	ibox64mod = 1 + int(rap / rs) + int(rbp / rs) * 4 + int(rcp / rs) * 4 * 4

	return
	end function
	

	integer*4 function ibox1000(raq, rbq, rcq)

	implicit none
	include "mybox.cmns"
	real*8 raq, rbq, rcq
	real*8 rap, rbp, rcp

	rap = raq - rcen2xtmp + 5.0D0 * rs
	rbp = rbq - rcen2ytmp + 5.0D0 * rs
	rcp = rcq - rcen2ztmp + 5.0D0 * rs

	ibox1000 = 1 + int(rap / rs) + int(rbp / rs) * 10 + int(rcp / rs) * 10 * 10

	return
	end function
