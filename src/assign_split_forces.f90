	subroutine assign_split_forces(iaq)

	implicit none
	include "mybox.cmns"

	integer*4 iaq, i, mn, m
	integer*4 jtype

	do i = 1, jsplitw2_total - 2, 3

		if(iw2force(i).eq.iaq) then

			mn = iw2force(i + 1)
			m  = iw2force(i + 2)

			!if(iaq.eq.12) write(*,'(f12.3,f12.3,f12.3)') fxw2(m, mn), fyw2(m, mn), fzw2(m, mn)
			!if(fxw2(m,mn).ne.0.0D0.or.fyw2(m,mn).ne.0.0D0.or.fzw2(m,mn).ne.0.0D0) write(*,*) "assign notzero"

			fxw2(m, mn) = fxw2(m, mn) + rw2force(i)
			fyw2(m, mn) = fyw2(m, mn) + rw2force(i + 1)
			fzw2(m, mn) = fzw2(m, mn) + rw2force(i + 2)

			jtype = m - 3 * ((m - 1) / 3)
			if(jtype.ne.2.and.jtype.ne.3) write(*,*) "assign_split_forces jtype"

			!if(iaq.eq.12) write(*,'(i5,i5,i5,i5)') iaq, iw2force(i), m, mn
			!if(iaq.eq.12) write(*,'(f12.3,f12.3,f12.3)') fxw2(m, mn), fyw2(m, mn), fzw2(m, mn)

		endif

	enddo

	do i = 1, jsplits_total - 2, 3

		if(isforce(i).eq.iaq) then

			mn = isforce(i + 1)
			m  = isforce(i + 2)

			!if(iaq.eq.12) write(*,'(f12.4,f12.4,f12.4)') fxs(m, mn), fys(m, mn), fzs(m, mn)
			!if(fxs(m,mn).ne.0.0D0.or.fys(m,mn).ne.0.0D0.or.fzs(m,mn).ne.0.0D0) write(*,*) "assign notzero"

			fxs(m, mn) = fxs(m, mn) + rsforce(i)
			fys(m, mn) = fys(m, mn) + rsforce(i + 1)
			fzs(m, mn) = fzs(m, mn) + rsforce(i + 2)

			!if(iaq.eq.12) write(*,'(i5,i5,i5,i5)') iaq, isforce(i), m, mn
			!if(iaq.eq.12) write(*,'(f12.4,f12.4,f12.4)') fxs(m, mn), fys(m, mn), fzs(m, mn)

		endif

	enddo

	!write(*,*) "split:", iaq, jsplitw2_total/3, jsplits_total/3

	return
	end subroutine
