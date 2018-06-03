using IterativeSolvers
include("usymqr.jl")
using Base.Test
using LinearMaps

function unsymmetric_problem(T, n)
    A = rand(T, n, n)
    x = ones(T, n)
    b = A * x
    A, x, b
end

srand(123)
n = 10

@testset "Unsymmetric Matrix{$T}" for T in (Float32, Float64)
    A, x, b = unsymmetric_problem(T, n)
    tol = sqrt(eps(real(T)))
    x0 = rand(T, n)

    x1, hist1 = usymqr(A, b, maxiter = 10n, tol = tol, log = true)

    @test isa(hist1, ConvergenceHistory)
    @test norm(b - A * x1) / norm(b) ≤ tol
    @test hist1.isconverged
end

@testset "SparseMatrixCSC{$T}" for T in (Float32, Float64)
    A = let
        B = sprand(n, n, 2 / n)
    end

    x = ones(T, n)
    b = A * x
    tol = sqrt(eps(real(T)))

    x_approx, hist = usymqr(A, b, maxiter = 10n, tol = tol, log = true)

    @test norm(b - A * x_approx) / norm(b) ≤ tol
    @test hist.isconverged
end
