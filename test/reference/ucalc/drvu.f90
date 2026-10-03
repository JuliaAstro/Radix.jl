!     Driver around XSTAR's real ucalc (see README.md). Reads cases from stdin:
!       ndesc nrdesc nrdt nidt            data type, rate type, number of reals / integers
!       (nidt integers)                    the record's integers (last one = ion index)
!       (nrdt reals)                       the record's reals
!       t xpx xee xh0 xh1 vturb trad cfrac abund1 abund2 ptmp1 ptmp2 amass
!                                          t in 1e4 K; densities in cm^-3; amass of the element
!       e0 dlnE index norm kpar gpar epar lfast radiation: epi(i)=e0 exp(dlnE (i-1)), bremsa(i)=norm epi^-index
!                                          (norm=0: no radiation); kpar, gpar: level index and statistical
!                                          weight and energy of the level of the next ion that photoionization leaves;
!                                          lfast is ucalc's speed switch (recombination integrals need >= 2)
!       nlev nspec                         levels of the ion, and how many are given below
!       (nspec lines) i E g Einf n 2S+1 L  the levels the record refers to
!     and prints ans1..ans6, idest1..idest4, opakab and four checksums of the
!     arrays ucalc fills (sums of opakc, opakcont, rccemis(1,:), rccemis(2,:)) for each case.
      program drvu
      use globaldata
      implicit none
      TYPE :: level_temp
        sequence
        real(8) :: rlev(10,ndl)
        integer:: ilev(10,ndl),nlpt(ndl),iltp(ndl)
        character(1) :: klev(100,ndl)
      END TYPE level_temp
      TYPE(level_temp) :: leveltemp
      integer ndesc,nrdesc,ml,lcon,jkion,indonly,nrdt,np1r,nidt,np1i
      integer nkdt,np1k,idest1,idest2,idest3,idest4,lpriu,ncn2,nlev
      integer lfast,lun11,i,ios,k,nspec
      real(8) vturbi,cfrac,ans1,ans2,ans3,ans4,ans5,ans6,abund1,abund2
      real(8) ptmp1,ptmp2,xpx,opakab,rr,delr,t,trad,tsq,xee,xh1,xh0
      real(8) E,g,einf,amass,e0,dlne,sindex,snorm,gpar,epar
      integer kpar
      integer n,l,s2
      character(49) kdesc2
      real(8), allocatable :: epi(:),bremsa(:),bremsint(:)
      real(8), allocatable :: rccemis(:,:),opakc(:),opakcont(:)
      allocate(epi(ncn),bremsa(ncn),bremsint(ncn),opakc(ncn),opakcont(ncn))
      allocate(rccemis(2,ncn))
      epi=0.d0; bremsa=0.d0; bremsint=0.d0; opakc=0.d0; opakcont=0.d0
      rccemis=0.d0
      allocate(derivedpointers%npcon(10),derivedpointers%npconi(10))
      allocate(derivedpointers%npilev(10,10),derivedpointers%npilevi(10))
      allocate(derivedpointers%npconi2(10))
      derivedpointers%npcon=0; derivedpointers%npconi=0
      derivedpointers%npilev=0; derivedpointers%npilevi=0
      derivedpointers%npconi2=0
      allocate(derivedpointers%npar(10),derivedpointers%nplini(10))
      allocate(derivedpointers%nplin(10),derivedpointers%npfi(60,10))
      allocate(derivedpointers%npfirst(60),derivedpointers%npnxt(10))
      derivedpointers%npar=0; derivedpointers%nplini=0
      derivedpointers%nplin=0; derivedpointers%npfi=0
      derivedpointers%npfirst=0; derivedpointers%npnxt=0
      if (.not.allocated(masterdata%nptrs)) allocate(masterdata%nptrs(10,10))
      masterdata%nptrs=0
      do i=1,ncn
        epi(i)=0.1d0*exp(0.0015d0*(i-1))
      enddo
      lun11=6; lcon=0; jkion=1; indonly=0; lfast=1; lpriu=0; ncn2=9999
      rr=1.d18; delr=1.d17; kdesc2=' '; ml=3; np1r=3; np1i=1; np1k=1; nkdt=0
 10   read (5,*,iostat=ios) ndesc,nrdesc,nrdt,nidt
      if (ios.ne.0) stop
      if (allocated(masterdata%idat1)) deallocate(masterdata%idat1)
      if (allocated(masterdata%rdat1)) deallocate(masterdata%rdat1)
      if (allocated(masterdata%kdat1)) deallocate(masterdata%kdat1)
      allocate(masterdata%idat1(nidt+20),masterdata%rdat1(nrdt+20),masterdata%kdat1(10))
      masterdata%idat1=0; masterdata%rdat1=0.d0
      read (5,*) (masterdata%idat1(i),i=1,nidt)
      read (5,*) (masterdata%rdat1(2+i),i=1,nrdt)
      read (5,*) t,xpx,xee,xh0,xh1,vturbi,trad,cfrac,abund1,abund2,ptmp1,ptmp2,amass
      ! pointer chain used by the line rates: record 3 (ml) -> ion record 2 ->
      ! element record 1, whose second real is the atomic mass
      masterdata%rdat1(1)=1.d0
      masterdata%rdat1(2)=amass
      masterdata%nptrs(1,1)=1; masterdata%nptrs(8,1)=1; masterdata%nptrs(9,1)=1
      masterdata%nptrs(10,1)=1
      masterdata%nptrs(1,3)=1; masterdata%nptrs(2,3)=ndesc
      masterdata%nptrs(3,3)=nrdesc; masterdata%nptrs(5,3)=nrdt
      masterdata%nptrs(6,3)=nidt; masterdata%nptrs(8,3)=3
      masterdata%nptrs(9,3)=1; masterdata%nptrs(10,3)=1
      derivedpointers%npar(3)=2; derivedpointers%npar(2)=1
      derivedpointers%nplini(3)=1; derivedpointers%nplin(1)=3
      tsq=sqrt(t)
      read (5,*) e0,dlne,sindex,snorm,kpar,gpar,epar,lfast
      do i=1,ncn
        epi(i)=e0*exp(dlne*(i-1))
        bremsa(i)=snorm*epi(i)**(-sindex)
      enddo
      opakc=0.d0; opakcont=0.d0; rccemis=0.d0; opakab=0.d0
      ! level record of the next ion (record 5) that the photoionization record leaves
      masterdata%idat1(nidt+1:nidt+6)=(/1,2,0,26,kpar,2/)
      masterdata%rdat1(2+nrdt+1:2+nrdt+4)=(/epar,gpar,1.d0,0.d0/)
      masterdata%nptrs(1,5)=1; masterdata%nptrs(2,5)=6; masterdata%nptrs(3,5)=13
      masterdata%nptrs(5,5)=4; masterdata%nptrs(6,5)=6; masterdata%nptrs(7,5)=0
      masterdata%nptrs(8,5)=2+nrdt+1; masterdata%nptrs(9,5)=nidt+1
      masterdata%nptrs(10,5)=1
      derivedpointers%npfi(13,2)=5; derivedpointers%npar(5)=4
      derivedpointers%npnxt(5)=6; derivedpointers%npar(6)=4   ! a real level always has a successor
      derivedpointers%npnxt(6)=0
      read (5,*) nlev,nspec
      do k=1,nspec
        read (5,*) i,E,g,einf,n,s2,l
        leveltemp%rlev(1,i)=E
        leveltemp%rlev(2,i)=g
        leveltemp%rlev(4,i)=einf
        leveltemp%ilev(1,i)=n
        leveltemp%ilev(2,i)=s2
        leveltemp%ilev(3,i)=l
      enddo
      call ucalc(ndesc,nrdesc,ml,lcon,jkion,vturbi,cfrac,indonly,       &
     &   nrdt,np1r,nidt,np1i,nkdt,np1k,ans1,ans2,                       &
     &   ans3,ans4,ans5,ans6,idest1,idest2,idest3,idest4,               &
     &   abund1,abund2,ptmp1,ptmp2,xpx,opakab,                          &
     &   opakc,opakcont,rccemis,lpriu,kdesc2,                           &
     &   rr,delr,t,trad,tsq,xee,xh1,xh0,                                &
     &   epi,ncn2,bremsa,bremsint,                                      &
     &   leveltemp,                                                     &
     &   nlev,lfast,lun11)
      write (6,'(6(1pe24.15),4i8,5(1pe24.15))') ans1,ans2,ans3,ans4,ans5,ans6,idest1,idest2,idest3,idest4,opakab, &
     &  sum(opakc),sum(opakcont),sum(rccemis(1,:)),sum(rccemis(2,:))
      goto 10
      end
