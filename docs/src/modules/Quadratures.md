# Quadratures

```@meta
CurrentModule = VortexPasta.Quadratures
CollapsedDocStrings = true
```

```@docs
Quadratures
```

## Quadrature rules

```@docs
AbstractQuadrature
```

### Fixed-size quadratures

These quadrature rules are meant to have a small size (typically less than 10).
They have virtually zero creation cost, i.e. doing `quad = GaussLegendre(4)` is basically free.

```@docs
StaticSizeQuadrature
GaussLegendre
NoQuadrature
```

### Variable-size quadratures

These quadrature rules should be constructed just once, as they allocate vectors.
These are usually adaptive quadratures.

```@docs
PreallocatedQuadrature
AdaptiveTanhSinh
```

## Estimating integrals

```@docs
Quadratures.integrate
```

## Computing quadrature rules

```@docs
quadrature
```

## Other functions

```@docs
length
```
