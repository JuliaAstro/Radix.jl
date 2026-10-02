#!/bin/sh
# Build a driver around XSTAR's real ucalc (Fortran 90 in ftools/xstar/xstarlib/src).
#
#   usage: build.sh /path/to/ftools/xstar/xstarlib/src  /path/to/builddir
#
# Needs gfortran (e.g. conda-forge). The flags mirror HEASoft's build (no
# -fdefault-real-8), so literals such as 0.861707 are single precision exactly as
# in the real xstar executable.
set -e
SRC=$1; OUT=$2; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT" && cd "$OUT"
FF="-c -O0 -w -ffree-line-length-none -fno-automatic -std=legacy -J$OUT -I$OUT -I$SRC"
for f in constants globaldata; do (cd "$SRC" && gfortran $FF $f.f90 -o "$OUT/$f.o"); done
for f in $(ls "$SRC"/*.f90 | xargs -n1 basename | sed 's/\.f90//'); do
  case $f in constants|globaldata|xstar|pprint) continue;; esac
  (cd "$SRC" && gfortran $FF $f.f90 -o "$OUT/$f.o")
done
gfortran -c -O0 -w -ffree-line-length-none -I"$OUT" -I"$SRC" -J"$OUT" "$HERE/drvu.f90" -o drvu.o
rm -f libx.a; ar rcs libx.a $(ls *.o | grep -v drvu.o)
gfortran -O0 -o drvu drvu.o libx.a
echo "built $OUT/drvu"
