# [Stationary methods](@id Stationary)

Stationary methods are typically used as smoothers in multigrid methods, where only very few iterations are applied to get rid of high-frequency components in the error. The implementations of stationary methods have this goal in mind, which means there is no other stopping criterion besides the maximum number of iterations.

!!! note "CSC versus CSR"
    Julia stores matrices column-major. In order to avoid cache misses, the implementations of our stationary methods traverse the matrices column-major. This deviates from classical textbook implementations. Also the SOR and SSOR methods cannot be computed efficiently in-place, but require a temporary vector.

    When it comes to `SparseMatrixCSC`, we precompute in all stationary methods an integer array of the indices of the diagonal to avoid expensive searches in each iteration.

## Jacobi

```@docs
jacobi
jacobi!
```

## Gauss-Seidel

```@docs
gauss_seidel
gauss_seidel!
```

## [Successive over-relaxation (SOR)](@id SOR)

```@docs
sor
sor!
```

## [Symmetric successive over-relaxation (SSOR)](@id SSOR)

```@docs
ssor
ssor!
```

!!! tip
    All stationary methods can be used a [iterators](@ref Iterators).