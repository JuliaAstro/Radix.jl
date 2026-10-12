Reference output from XSTAR 2.59j (HEASoft conda package, arm64) for the thick slab of `../xstar_thick_slab` with `npass=3` (the command of its README with `npass=3`; 10 minutes):
the first pass goes through the zones with the optical depths of the lines and edges from the source, the second comes back from the outer edge at the same temperatures to give the depths
beyond each zone, and the third goes through again with those as the outward depths of the escape probabilities.

- `xout_cont1.fits` and `xout_spect1.fits`: the incident and transmitted spectra and the emission columns (the same files as of the thick slab, after the three passes);
- `xout_lines1.fits` the table of the 599 strongest lines, and `xout_lines1_npass1.fits` the same of the run with `npass=1` (`../xstar_thick_slab`, whose `xout_lines1.fits` is not kept);
- the total of the table of the recombination edges (`xout_rrc1.fits`, 6 MB, not kept): 798877.05 (807150.61 with `npass=1`).

The lines change by up to a factor 16 (the 3p-3d line of Fe X) between the first pass and the third, and the He II 3s-2p line by 1.6: the escape probabilities of the thick lines (the He II Lyα has a depth of 830) change the
populations of their upper levels.
