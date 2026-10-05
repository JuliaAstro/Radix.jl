Reference output from XSTAR 2.59j (HEASoft conda package, arm64) for the first zone of a hot gas at log ξ = 1 whose iron is mostly
Fe VII-XII, with T = 4.5513 (10⁴ K) and x_e = 1 (the temperature of the equilibrium of `xeq` below, `niter=0` keeps `x_e`):

    xstar cfrac=0 temperature=4.5513 lcpres=0 pressure=0.03 density=1e4 spectrum=pow spectun=0 trad=-1 \
      rlrad38=1e8 column=1e17 rlogxi=1 abundtbl=xdef modelname=eq2 nsteps=1 niter=0 lwrite=0 \
      lprint=3 lstep=0 emult=1 taumax=5 xeemin=0.1 critf=1e-7 vturbi=0 radexp=0 ncn2=999 \
      loopcontrol=0 npass=1 mode=h

`xout_abund1.fits` (the ion fractions, first row) and `xout_cont1.fits` (the incident spectrum) are kept. The table of heating and cooling of `lprint=3`
(`xout_step.log`, 1.6 million lines) gives for iron 1.07571729E-15 and 9.67042360E-16 erg cm⁻³ s⁻¹ (photon point of view) and in total
1.18243531E-14 and 1.16261945E-14; the hydrogen of the database of the tree has a different cooling (see `../xstar_pow_xi2_eq`).

The same run with `niter=99`, `column=1e21` and `rlrad38=1e8` (the thick slab of `../xstar_thick_slab` with the iteration of the temperature, `emult=1`, `nsteps=3`) gives
T = 4.5513, 4.5513, 4.5513, 4.5513, 4.50231, 4.47051, 4.43977, 4.38722, 4.33672, 4.29065 (10⁴ K) in the rows of its table, `x_e` = 1.20772 and zones whose
inner edges are at 0, 0, 8.55289, 17.1058, 25.5628, 33.9275, 42.1408, 56.9934 and 71.5859 ×10¹⁵ cm. With the database of the XSTAR package (`$HEADAS/refdata/atdb.fits`) Radix
finds T = 4.5537 for the first zone and 4.2905 for the last (0.05% and 0.004%) and the same `x_e` and zones to 0.03%.
