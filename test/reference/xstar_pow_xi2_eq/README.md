Reference output from XSTAR 2.59j (HEASoft conda package, arm64) with the iteration of the temperature and the electron fraction
on: the same run as `../xstar_pow_xi2` but with `niter=99`.

    xstar cfrac=0 temperature=100 lcpres=0 pressure=0.03 density=1e4 spectrum=pow \
      spectun=0 trad=-1 rlrad38=1e-6 column=1e17 rlogxi=2 abundtbl=xdef \
      modelname=ref2 nsteps=2 niter=99 lwrite=0 lprint=0 lstep=0 emult=0.5 \
      taumax=5 xeemin=0.1 critf=1e-7 vturbi=1 radexp=0 ncn2=999 loopcontrol=0 \
      npass=1 mode=h

Only `xout_abund1.fits` is kept. The first zone (r = 10¹³ cm, log ξ = 2) settles at T = 18.8925 (10⁴ K) with `x_e` = 1.21001; the
radiation of the incident spectrum of `../xstar_pow_xi2/xout_cont1.fits` is the same.

The totals of the per-element table that XSTAR prints at `lprint=3` for the first zone at T = 10⁶ K and `x_e` = 1 (`../xstar_pow_xi2`
with `nsteps=1 lprint=3`): heating 5.20663621E-15, cooling 5.77911452E-15 erg cm⁻³ s⁻¹ (relative imbalance -0.104221974), Compton
2.37589869E-16 and 1.95382798E-16, free-free heating 2.44377767E-26, bremsstrahlung cooling 1.99246730E-16. Note that the database of
the XSTAR 2.59j package (`$HEADAS/refdata/atdb.fits`, 874.7 MB) is not the one of the XSTAR source tree (`ftools/xstar/data/atdb.fits`,
871.1 MB): for hydrogen the records of type 60 are not in it, and the 2s-3s record of the older database connects the levels 1 and 7.
