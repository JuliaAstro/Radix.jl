!     Driver around XSTAR's real linopac (the opacity of a line put into the continuum bins) and voigte.
!       stdin: ncn2, then the ncn2 energies epi (eV);
!              then nlines, and for each line: optpp, rcem1, rcem2, elin, vturbi, t, aatmp, delea, lfast
!       stdout: for each line the bins whose opacity or emissivity changed: index, opakc, rccemis(1), rccemis(2)
!       (opakc and rccemis start at zero for every line); then a table of voigte(v, a) for the pairs after the lines:
!       n ; then n pairs v, a.
!     Built like drvpesc (see build.sh and README.md): gfortran -c -O0 -w -ffree-line-length-none -fno-automatic -std=legacy ...
      program drvlinopac
      implicit none
      integer ncn
      parameter (ncn=999999)
      real(8), allocatable :: epi(:), opakc(:), rccemis(:,:)
      real(8) optpp, rcem1, rcem2, elin, vturbi, t, aatmp, delea, v, a, voigte
      integer ncn2, nlin, il, lfast, i, n, lprie, lun11
      allocate(epi(ncn), opakc(ncn), rccemis(2,ncn))
      epi = 0.d0
      read (5,*) ncn2
      read (5,*) (epi(i), i=1,ncn2)
      read (5,*) nlin
      lprie = 0
      lun11 = 6
      do il=1,nlin
        read (5,*) optpp, rcem1, rcem2, elin, vturbi, t, aatmp, delea, lfast
        opakc = 0.d0
        rccemis = 0.d0
        call linopac(lprie,lun11,optpp,rcem1,rcem2,elin,vturbi,t,aatmp,delea,epi,ncn2,opakc,rccemis,lfast)
        write (6,'(a,i6)') 'line', il
        do i=1,ncn
          if (opakc(i).ne.0.d0 .or. rccemis(1,i).ne.0.d0 .or. rccemis(2,i).ne.0.d0) &
     &      write (6,'(i8,3(1pe25.16))') i, opakc(i), rccemis(1,i), rccemis(2,i)
        enddo
      enddo
      read (5,*) n
      do i=1,n
        read (5,*) v, a
        write (6,'(a,3(1pe25.16))') 'voigt', v, a, voigte(v,a)
      enddo
      end
