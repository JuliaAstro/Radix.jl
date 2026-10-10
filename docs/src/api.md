# Reference

## Models and results

```@docs
Model
slab
march
tables
ConstantPressure
composition
```

## Sources

```@docs
powerlaw
blackbody
thermal
tabulated
```

## Atomic data

```@docs
Atom
AtomicLevel
AtomicLine
Levels
levels
counts
nlevels
weight
load
Radiation
Opacity
```

## Constants

```@docs
Constants
constants
constants!
ucalc
```

## Internals

The functions that are not exported (the rates of each record type, the level balance, the transfer and the continuum) have docstrings;
they are listed here and are called as `Radix.name`.

```@autodocs
Modules = [Radix]
Public = false
Order = [:type, :function]
```
