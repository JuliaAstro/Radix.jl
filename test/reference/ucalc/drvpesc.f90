!     Driver around XSTAR's real escape probabilities pescl (lines) and pescv (recombination continua).
!       stdin: n ; then n optical depths.  prints tau, pescl(tau), pescv(tau)
!     Built like drvu (see build.sh): gfortran -c -O0 -w -ffree-line-length-none -fno-automatic -std=legacy -I$OUT -J$OUT drvpesc.f90
!                                     gfortran -O0 -o drvpesc drvpesc.o $OUT/libx.a
      program drvpesc
      implicit none
      integer n,i
      real(8) tau,pescl,pescv
      read (5,*) n
      do i=1,n
        read (5,*) tau
        write (6,'(3(1pe25.16))') tau,pescl(tau),pescv(tau)
      enddo
      end
