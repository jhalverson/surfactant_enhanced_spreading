	subroutine evalw1w1(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq
	integer*4 i, i1, i2, j, j1, j2
	integer*4 ib4, ib216, ib64to216, ictk
	real*8 xOi, yOi, zOi, xOj, yOj, zOj
	real*8 xGi, yGi, zGi, xGj, yGj, zGj
	real*8 xHi, yHi, zHi, xHj, yHj, zHj
	real*8 xOOij, yOOij, zOOij, rOOijsq, rOOij
	real*8 xOGij, yOGij, zOGij, rOGijsq, rOGij
	real*8 xOHij, yOHij, zOHij, rOHijsq, rOHij
	real*8 xGOij, yGOij, zGOij, rGOijsq, rGOij
	real*8 xGGij, yGGij, zGGij, rGGijsq, rGGij
	real*8 xGHij, yGHij, zGHij, rGHijsq, rGHij
	real*8 xHOij, yHOij, zHOij, rHOijsq, rHOij
	real*8 xHGij, yHGij, zHGij, rHGijsq, rHGij
	real*8 xHHij, yHHij, zHHij, rHHijsq, rHHij
	real*8 fel, qOqO, qOqH, qHqH
	real*8 flj, ft, rsi, r6, rpl
	
	qOqO = rqw(1) * rqw(1) * xpi
	qOqH = rqw(1) * rqw(2) * xpi
	qHqH = rqw(2) * rqw(3) * xpi

	do ib4 = 1, 64

		ib216 = ib64to216(ib4)
		ictk = 3 * ictw1(ib216)

		do i = 1, ictk - 5, 3

			xOi = rxw1(i, ib216)
			yOi = ryw1(i, ib216)
			zOi = rzw1(i, ib216)

			i1 = i + 1
			xGi = rxw1(i1, ib216)
			yGi = ryw1(i1, ib216)
			zGi = rzw1(i1, ib216)

			i2 = i + 2
			xHi = rxw1(i2, ib216)
			yHi = ryw1(i2, ib216)
			zHi = rzw1(i2, ib216)

			do j = i + 3, ictk - 2, 3

				xOj = rxw1(j, ib216)
				yOj = ryw1(j, ib216)
				zOj = rzw1(j, ib216)

				j1 = j + 1
				xGj = rxw1(j1, ib216)
				yGj = ryw1(j1, ib216)
				zGj = rzw1(j1, ib216)

				j2 = j + 2
				xHj = rxw1(j2, ib216)
				yHj = ryw1(j2, ib216)
				zHj = rzw1(j2, ib216)

				xOOij = xOi - xOj
				yOOij = yOi - yOj
				zOOij = zOi - zOj

				rOOijsq = xOOij**2 + yOOij**2 + zOOij**2
				rOOij = sqrt(rOOijsq)

				fel = qOqO / rOOij**3

				flj = 0.0D0
				if(rOOij.le.rc) then
				
					rsi = 1.0D0 / rOOijsq
					r6 = rsi**3
					rpl = 48.0D0 * r6 * (r6 - 0.5D0)
					flj = rpl * rsi

				endif
				
				ft = fel + flj

				fxw1(i, ib216) = fxw1(i, ib216) + xOOij * ft
				fyw1(i, ib216) = fyw1(i, ib216) + yOOij * ft
				fzw1(i, ib216) = fzw1(i, ib216) + zOOij * ft

				fxw1(j, ib216) = fxw1(j, ib216) - xOOij * ft
				fyw1(j, ib216) = fyw1(j, ib216) - yOOij * ft
				fzw1(j, ib216) = fzw1(j, ib216) - zOOij * ft

				! do O with molecule j
				! O with G
				xOGij = xOi - xGj
				yOGij = yOi - yGj
				zOGij = zOi - zGj

				rOGijsq = xOGij**2 + yOGij**2 + zOGij**2
				rOGij = sqrt(rOGijsq)

				fel = qOqH / rOGij**3

				fxw1(i, ib216) = fxw1(i, ib216) + xOGij * fel
				fyw1(i, ib216) = fyw1(i, ib216) + yOGij * fel
				fzw1(i, ib216) = fzw1(i, ib216) + zOGij * fel

				fxw1(j1, ib216) = fxw1(j1, ib216) - xOGij * fel
				fyw1(j1, ib216) = fyw1(j1, ib216) - yOGij * fel
				fzw1(j1, ib216) = fzw1(j1, ib216) - zOGij * fel
				
				! O with H
				xOHij = xOi - xHj
				yOHij = yOi - yHj
				zOHij = zOi - zHj

				rOHijsq = xOHij**2 + yOHij**2 + zOHij**2
				rOHij = sqrt(rOHijsq)

				fel = qOqH / rOHij**3

				fxw1(i, ib216) = fxw1(i, ib216) + xOHij * fel
				fyw1(i, ib216) = fyw1(i, ib216) + yOHij * fel
				fzw1(i, ib216) = fzw1(i, ib216) + zOHij * fel

				fxw1(j2, ib216) = fxw1(j2, ib216) - xOHij * fel
				fyw1(j2, ib216) = fyw1(j2, ib216) - yOHij * fel
				fzw1(j2, ib216) = fzw1(j2, ib216) - zOHij * fel
				
				! do G with molecule j
				! G with O
				xGOij = xGi - xOj
				yGOij = yGi - yOj
				zGOij = zGi - zOj

				rGOijsq = xGOij**2 + yGOij**2 + zGOij**2
				rGOij = sqrt(rGOijsq)

				fel = qOqH / rGOij**3

				fxw1(i1, ib216) = fxw1(i1, ib216) + xGOij * fel
				fyw1(i1, ib216) = fyw1(i1, ib216) + yGOij * fel
				fzw1(i1, ib216) = fzw1(i1, ib216) + zGOij * fel

				fxw1(j, ib216) = fxw1(j, ib216) - xGOij * fel
				fyw1(j, ib216) = fyw1(j, ib216) - yGOij * fel
				fzw1(j, ib216) = fzw1(j, ib216) - zGOij * fel
				
				! G with G
				xGGij = xGi - xGj
				yGGij = yGi - yGj
				zGGij = zGi - zGj

				rGGijsq = xGGij**2 + yGGij**2 + zGGij**2
				rGGij = sqrt(rGGijsq)

				fel = qHqH / rGGij**3

				fxw1(i1, ib216) = fxw1(i1, ib216) + xGGij * fel
				fyw1(i1, ib216) = fyw1(i1, ib216) + yGGij * fel
				fzw1(i1, ib216) = fzw1(i1, ib216) + zGGij * fel

				fxw1(j1, ib216) = fxw1(j1, ib216) - xGGij * fel
				fyw1(j1, ib216) = fyw1(j1, ib216) - yGGij * fel
				fzw1(j1, ib216) = fzw1(j1, ib216) - zGGij * fel
				
				! G with H
				xGHij = xGi - xHj
				yGHij = yGi - yHj
				zGHij = zGi - zHj

				rGHijsq = xGHij**2 + yGHij**2 + zGHij**2
				rGHij = sqrt(rGHijsq)

				fel = qHqH / rGHij**3

				fxw1(i1, ib216) = fxw1(i1, ib216) + xGHij * fel
				fyw1(i1, ib216) = fyw1(i1, ib216) + yGHij * fel
				fzw1(i1, ib216) = fzw1(i1, ib216) + zGHij * fel

				fxw1(j2, ib216) = fxw1(j2, ib216) - xGHij * fel
				fyw1(j2, ib216) = fyw1(j2, ib216) - yGHij * fel
				fzw1(j2, ib216) = fzw1(j2, ib216) - zGHij * fel
				
				! do H with molecule j
				! H with O
				xHOij = xHi - xOj
				yHOij = yHi - yOj
				zHOij = zHi - zOj

				rHOijsq = xHOij**2 + yHOij**2 + zHOij**2
				rHOij = sqrt(rHOijsq)

				fel = qOqH / rHOij**3

				fxw1(i2, ib216) = fxw1(i2, ib216) + xHOij * fel
				fyw1(i2, ib216) = fyw1(i2, ib216) + yHOij * fel
				fzw1(i2, ib216) = fzw1(i2, ib216) + zHOij * fel

				fxw1(j, ib216) = fxw1(j, ib216) - xHOij * fel
				fyw1(j, ib216) = fyw1(j, ib216) - yHOij * fel
				fzw1(j, ib216) = fzw1(j, ib216) - zHOij * fel
				
				! H with G
				xHGij = xHi - xGj
				yHGij = yHi - yGj
				zHGij = zHi - zGj

				rHGijsq = xHGij**2 + yHGij**2 + zHGij**2
				rHGij = sqrt(rHGijsq)

				fel = qHqH / rHGij**3

				fxw1(i2, ib216) = fxw1(i2, ib216) + xHGij * fel
				fyw1(i2, ib216) = fyw1(i2, ib216) + yHGij * fel
				fzw1(i2, ib216) = fzw1(i2, ib216) + zHGij * fel

				fxw1(j1, ib216) = fxw1(j1, ib216) - xHGij * fel
				fyw1(j1, ib216) = fyw1(j1, ib216) - yHGij * fel
				fzw1(j1, ib216) = fzw1(j1, ib216) - zHGij * fel
				
				! H with H
				xHHij = xHi - xHj
				yHHij = yHi - yHj
				zHHij = zHi - zHj

				rHHijsq = xHHij**2 + yHHij**2 + zHHij**2
				rHHij = sqrt(rHHijsq)

				fel = qHqH / rHHij**3

				fxw1(i2, ib216) = fxw1(i2, ib216) + xHHij * fel
				fyw1(i2, ib216) = fyw1(i2, ib216) + yHHij * fel
				fzw1(i2, ib216) = fzw1(i2, ib216) + zHHij * fel

				fxw1(j2, ib216) = fxw1(j2, ib216) - xHHij * fel
				fyw1(j2, ib216) = fyw1(j2, ib216) - yHHij * fel
				fzw1(j2, ib216) = fzw1(j2, ib216) - zHHij * fel
				
			enddo
			
		enddo

	enddo

	return
	end subroutine
