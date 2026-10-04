!     Driver around XSTAR's real comp2 (and cmpfnc): the Compton heating and cooling integrals of a spectrum.
!       usage:  drvcomp /path/to/coheat.dat < input
!       input:  t n                 temperature (1e4 K) and number of energies
!               epi bremsa          n lines: energy (eV) and flux (erg/s/cm^2/erg)
!     prints cmp1 and cmp2 (heating = cmp1 nx ergsev, cooling = kT cmp2 nx ergsev). Built like drvu (see build.sh):
!       gfortran -c -O0 -w -ffree-line-length-none -fno-automatic -std=legacy -I$OUT -J$OUT drvcomp.f90
!       gfortran -O0 -o drvcomp drvcomp.o $OUT/libx.a
      program drvcomp
      use globaldata
      implicit none
      integer n,i,mm,ll,lpri,lun11
      real(8) t,cmp1,cmp2
      real(8), allocatable :: epi(:),bremsa(:)
      character(256) datafil3
      call get_command_argument(1,datafil3)
      open(unit=25,file=datafil3,status='old')
      do mm=1,ncomp
        do ll=1,ncomp
          read (25,902) sxcomp(mm),ecomp(ll),decomp(mm,ll)
          enddo
        enddo
      close(25)
  902 format (9x,e12.4,12x,2(e12.4))
      read (5,*) t,n
      allocate(epi(ncn),bremsa(ncn))
      epi=0.d0; bremsa=0.d0
      do i=1,n
        read (5,*) epi(i),bremsa(i)
      enddo
!     (comp2 and cmpfnc assign to lpri: it must be a variable)
      lpri=0; lun11=6
      call comp2(lpri,lun11,epi,n,bremsa,t,cmp1,cmp2)
      write (6,'(2(1pe25.16))') cmp1,cmp2
      end
