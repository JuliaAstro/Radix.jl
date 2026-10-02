# Fortran reference for ucalc type 63 (CollisionProb)

`../collisionprob_cases.txt` holds `ans1`, `ans2` computed by XSTAR's own Fortran for
800 level pairs sampled from the real database at four (T, nₑ) conditions. It
was built with the double-precision build of the extracted routines, so the
Julia port can be checked to ~1e-9. XSTAR's single-precision build differs from
these values by up to 0.7%.

To regenerate (needs a Fortran compiler, e.g. `gfortran`):

```bash
python3 extract_type63.py /path/to/xstarsub.f
gfortran -std=legacy -w -fno-automatic -fdefault-real-8 -fdefault-double-8 \
    -o ref63 drv.f prob63.f routines63.f
./ref63 < cases.txt      # lines: ni li nf lf Z T[K] E1 E2 g1 g2 ne[cm-3]
```
