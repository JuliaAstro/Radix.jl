Reference output from XSTAR 2.59j (HEASoft conda package, arm64) for the thick slab of `../xstar_thick_slab` with hydrogen, helium, carbon and argon only (the other abundances set to 0), with `npass=1` and
`npass=3` and `lwrite=1` (the commands of `../xstar_hhe` with `cabund` and `arabund` left at 1 and `npass=<n> lwrite=1`):

- `lines_npass1.csv`, `lines_npass3.csv`: the lines of C V and Ar IX of `xout_lines1.fits` (`emit_inward` and `emit_outward` in 10³⁸ erg/s);
- `populations_zone1.csv`, `populations_zone3.csv`: the populations of the levels of C and Ar in the first and the third zones (`xo01_detail.fits`, written in the first pass; the third pass has the same).

Between the first pass and the third XSTAR changes the forbidden line of C V at 41.47 Å (2³S) by +10.5% and the 3s ¹P line of Ar IX at 48.74 Å by +13.7% (Radix: -0.06% and +4.4%). Both upper levels have populations below 10⁻⁶ of the
element (2³S of C V 1.2×10⁻⁹ in the first zone, 2.1×10⁻¹⁰ in the third; the 3s levels of Ar IX 10⁻¹¹ and less), which the relaxation of `msolvelucy` leaves out, and XSTAR's populations are not a steady state
of its own matrix (the matrix of Radix at the zone): the 2³S of C V in the first zone gets 6% of what it loses (Radix's steady state is 7 times lower), and in the third zone 3% less than it loses (18% above Radix); the 3³P of Ar IX in the third zone
are 13-23% above Radix's (residual 3-7%), while in the first zone they agree to 0.2%. See `docs/src/field-names.md`.
