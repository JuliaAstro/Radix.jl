      program drv
      implicit none
      integer idat(10),ilev(6,3),ni,li,nf,lf,iz,i
      real rlev(4,3),t,xnx,ans1,ans2,e1,e2,g1,g2
      integer ios
 10   read (5,*,iostat=ios) ni,li,nf,lf,iz,t,e1,e2,g1,g2,xnx
      if (ios.ne.0) stop
      do i=1,10
        idat(i)=0
      enddo
      idat(1)=1
      idat(2)=1
      idat(3)=2
      idat(4)=iz
      idat(5)=7
      ilev(1,1)=ni
      ilev(3,1)=li
      ilev(1,2)=nf
      ilev(3,2)=lf
      rlev(1,1)=e1
      rlev(1,2)=e2
      rlev(2,1)=g1
      rlev(2,2)=g2
      call prob63(idat,5,ilev,rlev,2,t/1.e4,xnx,0,6,ans1,ans2)
      write (6,'(1p2e24.15)') ans1,ans2
      goto 10
      end
