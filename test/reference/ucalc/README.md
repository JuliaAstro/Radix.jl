# Reference values from XSTAR's real `ucalc`

`../../ucalc_tests.jl` compares the Julia rates with what the actual Fortran 90
`ucalc` (ftools/xstar/xstarlib/src) returns for records of `atdb.fits`.

- `drvu.f90` fills XSTAR's module data for one database record (integers, reals,
  level data, a small pointer chain to an element record) and calls `ucalc`.
- `build.sh` compiles the xstarlib sources with `gfortran` the way HEASoft does
  (single-precision literals included) and links the driver.
- `generate.jl` samples real records of one data type, with their ion's level data,
  and writes the driver input.
- `data/typeNN.in` and `data/typeNN.out` are the saved inputs and `ucalc` outputs
  (`ans1..ans6, idest1..idest4, opakab`); `type51n9` is a synthetic 9-point CHIANTI
  case.

To regenerate a type (needs `gfortran` and the XSTAR sources):

```bash
./build.sh /path/to/ftools/xstar/xstarlib/src /tmp/ucalc-build
julia generate.jl /path/to/atdb.fits /tmp/ucalc-ref 56          # writes type56.in
/tmp/ucalc-build/drvu < /tmp/ucalc-ref/type56.in > data/type56.out
```

`ucalc`'s `t` is in units of 10⁴ K. The input format is documented at the top of
`drvu.f90` and in `generate.jl`.

`drvcomp.f90` is a second driver, around `comp2` (the Compton heating and cooling integrals of the spectrum that
`src/heating.jl` ports as `integral(Compton, …)`): it reads `coheat.dat` as `xstarsetup` does, then a temperature and a spectrum. For
the incident spectrum of `../xstar_pow_xi2` at 10¹³ cm (`F = L/(4π r²)`, 999 points) and T = 10⁶ K it returns
`cmp1 = 9.8459456377190043E-09` and `cmp2 = 8.9007105172073563E-11`, which `test/mixture_tests.jl` compares with.

`drvpesc.f90` is a third driver, around `pescl` and `pescv` (the escape probabilities of the lines and the recombination continua that
`src/transfer.jl` ports): it prints both for the optical depths on its input; `test/transfer_tests.jl` has the values for 14 depths.

`drvlinopac.f90` is a fourth driver, around `linopac` (the opacity of a line put into the continuum bins) and `voigte` that
`src/lines.jl` ports as `add_line!` and `voigt`: it reads a grid, lines (`optpp`, `rcem1`, `rcem2`, `elin`, `vturbi`, `t`, `aatmp`, `delea`, `lfast`) and pairs `(v, a)`, and
prints the bins that each line changed and `voigte`. `data/linopac.in` (a log grid of 999 points and four lines, two of them with `lfast = 2` and a damping parameter
below and above 10⁻⁶, one with `lfast = 3` and a hydrogen line) and `data/linopac.out` are what `test/lines_tests.jl` compares with. Built like `drvpesc`:
`gfortran -c -O0 -w -ffree-line-length-none -fno-automatic -std=legacy -I$OUT -J$OUT drvlinopac.f90`, then `gfortran -O0 -o drvlinopac drvlinopac.o $OUT/libx.a`.
