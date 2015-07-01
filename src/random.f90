	subroutine random(x,xr)

	implicit real*8 (a-h,o-z)
	data k,j,m,rm/5701,3612,566927,566927./
	ix=int(x*rm)
	kr=j*ix+k
	irand=mod(kr,m)
	xr=(float(irand)+0.5D0)/rm
	return
	end subroutine