Reference output from XSTAR 2.59j (HEASoft conda package, arm64) with the density from a file (`radexp=-100`: `density.dat`, a line of `r n` for each zone;
40 lines of a falling density with steps every five zones, from 3.16×10¹⁴ to 4.4×10¹⁵ cm), a power law of 10³⁵ erg/s, `log ξ = 2` at the first line, T = 10⁶ K and x_e = 1 (`niter=0`), without turbulence:

    xstar cfrac=0 temperature=100 lcpres=0 pressure=0.03 density=1e4 spectrum=pow spectun=0 trad=-1 rlrad38=1e-3 column=<N> rlogxi=2 abundtbl=xdef \
      modelname=src nsteps=2 niter=0 lwrite=0 lprint=0 lstep=0 emult=1 taumax=5 xeemin=0.1 critf=1e-7 vturbi=0 radexp=-100 ncn2=999 loopcontrol=0 npass=1 mode=h

- `xout_abund1_column.fits`: `column=1e18`. The table has six rows: the zones at the first five lines, and a last row at the sixth line past them (4.7449×10¹⁴ cm, 4.62×10³ cm⁻³:
  the line where the column passes 10¹⁸ cm⁻²) with the ion fractions of the zone before (4.476×10¹⁴ cm, 4.16×10³ cm⁻³). The distance, density and ionization parameter of
  that row are those of the final calculation of XSTAR, which `lprint=3` prints: for it the heating of the photons is 1.14746613E-15 and the cooling 1.26292604E-15 erg cm⁻³ s⁻¹ (totals of its table).
- `xout_abund1_table.fits`: `column=1e19`, more than the table has (4.4×10¹⁸): the zones are the 40 lines and the run ends with the file.

The zones are at the distances of the file whatever the thickness that `step` gives (that one is what the zone has for its opacity and emission: the distance of the next zone is the next line) and the column
grows with the density of the next line times its distance from the zone.
