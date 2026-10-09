Reference output from XSTAR 2.59j (HEASoft conda package, arm64) for a slab that is thick enough to show the absorption of the
continuum: a plane-parallel slab (the source is at 3.2×10²⁰ cm and the slab is 10¹⁷ cm thick) of column density 10²¹ cm⁻², at
log ξ = 1, T = 10⁵ K and an electron fraction of 1, both kept (`niter=0`).

    xstar cfrac=0 temperature=10 lcpres=0 pressure=0.03 density=1e4 spectrum=pow spectun=0 trad=-1 \
      rlrad38=1e8 column=1e21 rlogxi=1 abundtbl=xdef modelname=thick2 nsteps=3 niter=0 lwrite=0 \
      lprint=0 lstep=0 emult=1 taumax=5 xeemin=0.1 critf=1e-7 vturbi=0 radexp=0 ncn2=999 \
      loopcontrol=0 npass=1 mode=h

`xout_cont1.fits` (the incident and transmitted spectra) and `xout_abund1.fits` (the ion fractions of the zones: the rows are the
states at the inner edges of the zones, `delta_r` is the depth of the zone edge) are kept.

**`vturbi=0` matters.** The same run with `vturbi=1` (the default, `xout_cont1_vturb1.fits`) has a transmitted spectrum that is exactly the
incident one between 0.15 eV and 20 keV (and Thomson scattering outside it), and diffuse emission of 10⁻¹⁷ instead of 10¹⁵. The cause is
`gsmooth2` (called when `vturbi > 1e-34`): its Gaussian kernel is `E vtherm/c` wide, which is much narrower than the 1.6% bins of a grid of
999 points (and than the 0.16% bins of 9999 for T below 4×10⁵ K), the loop ends at the first neighbour with `exp(-earg) = 0`
and never adds the bin itself, so the smoothed opacity and emissivity arrays are 0 below 20 keV. The continuum is then not absorbed by the
slab at all. Radix does not do it.

`xout_abund1_npass3.fits` is the same run with `npass=3` (it takes four times as long): the ion fractions of the zones
after the pass back from the outer edge, which `passes=3` of `slab_model` reproduces.
