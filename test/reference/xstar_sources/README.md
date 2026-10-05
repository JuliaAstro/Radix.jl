The incident spectra (the `incident` column of `xout_cont1.fits`, in 10³⁸ erg/s/erg) that XSTAR 2.59j (HEASoft conda package, arm64) makes for its three thermal spectra, for
`rlrad38=1e-3` (a luminosity of 10³⁵ erg/s) and 999 energies; `trad` is the temperature in 10⁷ K:

    xstar cfrac=0 temperature=100 lcpres=0 pressure=0.03 density=1e4 spectrum=<type> spectun=0 trad=<T> \
      rlrad38=1e-3 column=1e17 rlogxi=2 abundtbl=xdef modelname=src nsteps=1 niter=0 lwrite=0 \
      lprint=0 lstep=0 emult=1 taumax=5 xeemin=0.1 critf=1e-7 vturbi=0 radexp=0 ncn2=999 \
      loopcontrol=0 npass=1 mode=h

- `bbody_0.5.fits`: `spectrum=bbody trad=0.5` (kT = 431 eV);
- `bbody_0.0003.fits`: `spectrum=bbody trad=0.0003`, for which x = E/kT is above 150 in most of the grid (the tail 1/E of `starf`, with the luminosity 1 that `xstar` gives it);
- `brems_0.1.fits`: `spectrum=brems trad=0.1` (kT = 86.2 eV).
