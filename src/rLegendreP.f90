      ! This function was taken from www.nr.com
      ! The first two if statements are modifications.
      real*8 FUNCTION rLegendreP(l, m, xx)

      implicit none
      INTEGER l,m
      real*8 xx
      INTEGER i,ll
      real*8 fact,pll,pmm,pmmp1,somx2
      
      if(l.eq.-1.and.m.eq.0) then
         rLegendreP = 1.0D0
         return
      endif
      
      if(m.gt.l) then
         rLegendreP = 0.0D0
         return
      endif
      
      !if(m.lt.0.or.m.gt.l.or.abs(xx).gt.1.) pause 'bad arguments in rLegendreP'
      pmm=1.
      if(m.gt.0) then
         somx2=sqrt((1.-xx)*(1.+xx))
         fact=1.
         do i=1,m
            pmm=-pmm*fact*somx2
            fact=fact+2.
         enddo
      endif
      if(l.eq.m) then
         rLegendreP=pmm
      else
         pmmp1=xx*(2*m+1)*pmm
         if(l.eq.m+1) then
            rLegendreP=pmmp1
         else
            do ll=m+2,l
               pll=(xx*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
               pmm=pmmp1
               pmmp1=pll
            enddo
            rLegendreP=pll
         endif
      endif
      
      return
      END FUNCTION
