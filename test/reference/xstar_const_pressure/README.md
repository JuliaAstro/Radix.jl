Reference output from XSTAR 2.59j (HEASoft conda package, arm64) for a slab at constant pressure with the iteration of the temperature
(`xout_abund1.fits` is kept):

    xstar cfrac=0 temperature=10 lcpres=1 pressure=1e-7 density=1e4 spectrum=pow spectun=0 trad=-1 \
      rlrad38=1e8 column=1e21 rlogxi=1 abundtbl=xdef modelname=cp1 nsteps=3 niter=99 lwrite=0 \
      lprint=0 lstep=0 emult=1 taumax=5 xeemin=0.1 critf=1e-7 vturbi=0 radexp=0 ncn2=999 \
      loopcontrol=0 npass=1 mode=h

With `lcpres=1` the hydrogen density is `P/(1.38×10⁻¹² T)` (T in 10⁴ K) at the temperature of each zone, and `rlogxi` is the logarithm of the ionization parameter of the pressure
`Ξ = L/(4π c P r²)`, so that r = 1.62965×10²⁰ cm. The first zone settles at T = 10.4611, n = 6927 and x_e = 1.20938; the table has three rows
(the zones at 0, 0 and 1.44363×10¹⁷ cm). Radix with the database of the package finds 10.4574, 6929 and 1.20938 for the first zone (0.04%, 0.03%).

The second zone of XSTAR keeps T = 10.4611 with a heating imbalance of 0.15% (its table): the iteration of its `dsec` leaves the temperature where it is when the false
position repeats a value (`testt < 2×10⁻⁹`). Radix with the same state finds the imbalance 0.154% and its converged temperature for that zone is 10.71.
