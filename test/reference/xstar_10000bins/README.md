Reference output from XSTAR 2.59j (HEASoft conda package, arm64) with `ncn2=10000` energy bins (the references of the other directories have 999; XSTAR's default is 9999):

- `thin_abund1.fits`: the first zone of the run of `../xstar_pow_xi2_eq` (log ξ = 2, 10⁴ cm⁻³, 10³² erg/s, T and x_e iterated: 18.8982) with `column=1e17 nsteps=2 emult=0.5 niter=99 ncn2=10000 vturbi=0`;
- `thick_abund1.fits`, `thick_cont1.fits`: the thick slab of `../xstar_thick_slab` (10²¹ cm⁻², log ξ = 1, T = 10⁵ K, `niter=0 nsteps=3 emult=1 vturbi=0`) with `ncn2=10000`: 15 rows (the zones of XSTAR's `step` are not those of Radix after the eighth: the
  `taumax = 5` condition of its rule includes the bin of the O VII resonance line at 574.8 eV after the sixth zone with a depth of 4.9876, and Radix's depth there is 5.0252, so Radix leaves that bin out of the minimum; the
  opacities of the zone agree to 0.75% in that bin and to 10⁻⁵ in the median).

XSTAR solves the balance, the heating and the cooling, and the edges on a grid of 999 bins (`ncn2m`, mapped by `bremsmap`) whatever `ncn2` is; only the continuum arrays (`calc_emis_all`) are on the grid of `ncn2` bins. The first zone of the thin slab,
T = 18.8925 at 999 bins and 18.8982 at 10000, changes by 0.03%.
