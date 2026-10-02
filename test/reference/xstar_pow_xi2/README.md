Reference output from XSTAR 2.59j (HEASoft conda package, arm64).

Command (run in an empty directory with PFILES pointing at it):

    xstar cfrac=0 temperature=100 lcpres=0 pressure=0.03 density=1e4 spectrum=pow \
      spectun=0 trad=-1 rlrad38=1e-6 column=1e17 rlogxi=2 abundtbl=xdef \
      modelname=ref1 nsteps=2 niter=0 lwrite=1 lprint=0 lstep=0 emult=0.5 \
      taumax=5 xeemin=0.1 critf=1e-7 vturbi=1 radexp=0 ncn2=999 loopcontrol=0 \
      npass=1 mode=h

The large xo01_detail*.fits files (up to 137 MB) are not kept.
