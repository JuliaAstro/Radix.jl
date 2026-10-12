Reference output from XSTAR 2.59j (HEASoft conda package, arm64) for the thick slab of `../xstar_thick_slab` with hydrogen and helium only (the abundances of
the other 28 elements set to 0), so that its tables of the lines list every line (with the metals the table has the strongest 600):

    xstar cfrac=0 temperature=10 lcpres=0 pressure=0.03 density=1e4 spectrum=pow spectun=0 trad=-1 rlrad38=1e8 column=1e21 rlogxi=1 abundtbl=xdef \
      modelname=hhe nsteps=3 niter=0 lwrite=0 lprint=0 lstep=0 emult=1 taumax=5 xeemin=0.1 critf=1e-7 vturbi=0 radexp=0 ncn2=999 loopcontrol=0 npass=1 mode=h \
      liabund=0 beabund=0 babund=0 cabund=0 ... znabund=0

- `lines.csv`: the table of the lines of `xout_lines1.fits` (435 lines, `emit_inward` and `emit_outward` in 10³⁸ erg/s);
- `populations_zone1.csv`: the populations of the levels of H I, He I and He II in the first zone (`lwrite=1`, `xo01_detail.fits`, relative to the element: the continuum is
  the next ion; the table has 112 levels).

The He I lines of the triplets (10830 Å, 3889, 5876, 7065, ...) are 2 times weaker in XSTAR than in Radix, and the population of 2³S is 5.4 times lower (9.75×10⁻¹¹ against 5.25×10⁻¹⁰). The
matrix of `calc_hmc_element` at `lprint=5` is the same as Radix's to 10⁻³ (the loss of 2³S is 0.0435213 and 0.0435308 s⁻¹), but XSTAR's populations are **not a steady state of it**: the flow into 2³S is 1.53 times its flow
out (the inflow is 9.05×10⁻¹² and the outflow 4.24×10⁻¹²), and one Jacobi sweep of that matrix from XSTAR's populations raises 2³S by 2.1. `msolvelucy` relaxes only
the populations above `eps = 10⁻⁶` and the levels of He I are all below it (the ground state is 5×10⁻⁸ of the helium). See `docs/src/field-names.md`.
