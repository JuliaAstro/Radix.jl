Reference output from XSTAR 2.59j (HEASoft conda package, arm64, database of the package) for a covering fraction of 0.5 with the temperature iterated (`niter=99`):

    xstar cfrac=0.5 temperature=100 lcpres=0 pressure=0.03 density=1e4 spectrum=pow spectun=0 trad=-1 \
      rlrad38=1e-2 column=1e19 rlogxi=2 abundtbl=xdef modelname=v nsteps=3 niter=99 lwrite=0 lprint=0 \
      lstep=0 emult=1 taumax=5 xeemin=0.1 critf=1e-7 vturbi=0 radexp=0 ncn2=999 loopcontrol=0 npass=1 mode=h

Only `xout_abund1.fits` is kept: four rows (the first zone has no thickness and appears twice; the last is the end of the slab) at r = 10¹⁵ cm. The covering fraction reduces the
photoexcitation of the lines (`ucalc` type 50) and the escape probabilities and Thomson scattering: the first zone is at T = 18.081 (10⁴ K), against 18.89 without it
(`../xstar_pow_xi2_eq`, the same ionization parameter and density at another luminosity and distance).
