# Intermolecular internuclear coupling data structures and methods

## Data structures

```@docs
DeltaKIMolBlockData
DeltaKISystemData
```

## Morse potential harmonic coupling function

This function can be used independently of AARMBEM-C.jl.

```@docs
keffMorse
```

## Functions to compute ``\Delta K_{\mathrm{I}}``

```@docs
readDeltaKISystemData
assembleDeltaKI!
```